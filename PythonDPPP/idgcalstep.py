from dppp import DPStep

import time
import idg

import idg.util as util

import astropy.io.fits as fits
import scipy.linalg
from .idgwindow import idgwindow
import numpy as np
import h5py

def next_composite(n):
    n += (n & 1)
    while True:
        nn = n
        while ((nn % 2) == 0): nn /= 2
        while ((nn % 3) == 0): nn /= 3
        while ((nn % 5) == 0): nn /= 5
        if (nn == 1): return n
        n += 2

class IDGCalStep(DPStep):

    def __init__(self,parset, prefix):
        super().__init__()
        self.parset = parset
        self.dpbuffers = []
        print(dir(parset))
        print("prefix = ", prefix)
        print(parset.getString(prefix + "type"))

        self.model_image = parset.getString(prefix + "modelimage")
        self.phase_interval = parset.getInt(prefix + "phaseinterval")
        self.amplitude_interval = parset.getInt(prefix + "amplitudeinterval")

        self.nr_timeslots = (self.amplitude_interval + self.phase_interval - 1) // self.phase_interval
        self.nr_timesteps_per_slot = self.phase_interval
        self.nr_timesteps = self.nr_timeslots * self.nr_timesteps_per_slot

        self.padding = parset.getFloat(prefix + "padding")
        self.nr_parameters_ampl = parset.getInt(prefix + "nrparametersamplitude")
        self.nr_parameters_phase = parset.getInt(prefix + "nrparametersphase")
        self.nr_parameters0 = self.nr_parameters_ampl + self.nr_parameters_phase
        self.nr_parameters = self.nr_parameters_ampl + self.nr_parameters_phase*self.nr_timeslots

        self.idg_parm = parset.getString(prefix + "idgparm")
        self.subgrid_size = 32
        self.nr_correlations = 4

        self.taper_support    = 7
        self.wterm_support    = 5
        self.aterm_support    = 5
        self.kernel_size      = self.taper_support + self.wterm_support + self.aterm_support

        self.solver_update_gain = 0.5
        self.pinv_tol = 1e-9

        self.initialized = False

    def show(self) :
        print("MyDPStep.show")

    def process(self, dpbuffer) :

        if not self.initialized:
            self.init()

        self.dpbuffers.append(dpbuffer)
        if len(self.dpbuffers) == self.amplitude_interval:
            self.process_buffers()
            self.dpbuffers = []


    def finish(self):
        self.process_buffers()
        self.output_file.close()


    def process_buffers(self):
        self.calibrate()
        for dpbuffer in self.dpbuffers:
            pass
            self.process_next_step(dpbuffer)

    def init(self):
        self.initialized = True

        self.proxy = idg.HybridCUDA.GenericOptimized(self.nr_correlations, self.subgrid_size)

        self.init_buffers()
        self.init_grid()
        self.init_idg_h5parm()
        self.init_basefunctions()


    def fill_buffers(self):
        for i, dpbuffer in enumerate(self.dpbuffers):
            self.uvw[:, i]['u'] = np.array(dpbuffer.get_uvw())[:,0]
            self.uvw[:, i]['v'] = -np.array(dpbuffer.get_uvw())[:,1]
            self.uvw[:, i]['w'] = -np.array(dpbuffer.get_uvw())[:,2]
            self.visibilities[:,i,:,:] = dpbuffer.get_data()
            self.weights[:,i,:,:] = dpbuffer.get_weights()
            # TODO
            # flags
            # vis_block[np.isnan(vis_block)] = 0

    def calibrate(self):


        self.fill_buffers()

        # Grid visibilities
        w_step = 400.0

        shift = np.array((0.0, 0.0, 0.0), dtype=np.float32)

        aterms         = util.get_identity_aterms(
                            self.nr_timeslots, self.nr_stations, self.subgrid_size, self.nr_correlations)
        aterms_offsets = util.get_example_aterms_offset(
                            self.nr_timeslots, self.nr_timesteps)



        print("calibrate_init...")
        print(self.visibilities.shape)
        print(self.weights.shape)
        print(self.uvw.shape)
        self.proxy.calibrate_init(
            w_step,
            shift,
            self.cell_size,
            self.kernel_size,
            self.subgrid_size,
            self.frequencies,
            self.visibilities,
            self.weights,
            self.uvw,
            self.baselines,
            self.grid,
            aterms_offsets,
            self.taper2)
        print("done.")

        X0 =  np.ones((self.nr_stations, 1))
        X1 =  np.zeros((self.nr_stations, 1))

        parameters = np.concatenate((X0,) + (self.nr_parameters-1)*(X1,), axis=1)

        aa = self.Bampl.sum(axis=1).sum(axis=1)
        parameters[:,:] = np.concatenate((aa[:self.nr_parameters_ampl,0], np.zeros((self.nr_parameters_phase*self.nr_timeslots,))))[np.newaxis, :]

        aterm_ampl = np.tensordot(parameters[:,:self.nr_parameters_ampl] , self.Bampl, axes = ((1,), (0,)))
        aterm_phase = np.exp(1j*np.tensordot(parameters[:,self.nr_parameters_ampl:].reshape((self.nr_stations, self.nr_timeslots, self.nr_parameters_phase)), self.Bphase, axes = ((2,), (0,))))
        aterms[:,:,:,:,:] = aterm_phase.transpose((1,0,2,3,4))*aterm_ampl

        nr_iterations = 0
        converged = False

        max_dx = 0.0

        timer = -time.time()

        timer0 = 0
        timer1 = 0

        previous_residual = 0.0

        while True:

            nr_iterations += 1

            #if nr_iterations == 2:
                #np.save("p2", parameters)
                #break

            print("iteration nr {0} ".format(nr_iterations),)

            max_dx = 0.0
            norm_dx = 0.0
            residual_sum = 0.0
            for i in range(self.nr_stations):
                timer1 -= time.time()

                # predict visibilities for current solution

                hessian  = np.zeros((self.nr_timeslots, self.nr_parameters0, self.nr_parameters0), dtype = np.float64)
                gradient = np.zeros((self.nr_timeslots, self.nr_parameters0), dtype = np.float64)
                residual = np.zeros((1, ), dtype = np.float64)

                aterm_ampl = np.repeat(np.tensordot(parameters[i,:self.nr_parameters_ampl], self.Bampl, axes = ((0,), (0,)))[np.newaxis,:], self.nr_timeslots, axis=0)
                aterm_phase = np.exp(1j * np.tensordot(parameters[i,self.nr_parameters_ampl:].reshape((self.nr_timeslots, self.nr_parameters_phase)), self.Bphase, axes = ((1,), (0,))))

                aterm_derivatives_ampl = aterm_phase[:, np.newaxis,:,:,:]*self.Bampl[np.newaxis,:,:,:,:]

                aterm_derivatives_phase = 1j*aterm_ampl[:, np.newaxis,:,:,:] * aterm_phase[:, np.newaxis,:,:,:] * self.Bphase[np.newaxis,:,:,:,:]

                aterm_derivatives = np.concatenate((aterm_derivatives_ampl, aterm_derivatives_phase), axis=1)
                aterm_derivatives = np.ascontiguousarray(aterm_derivatives, dtype=np.complex64)

                timer0 -= time.time()
                self.proxy.calibrate_update(i, aterms, aterm_derivatives, hessian, gradient, residual)

                timer0 += time.time()

                residual0 = residual[0]

                residual_sum += residual[0]

                gradient = np.concatenate((np.sum(gradient[:,:self.nr_parameters_ampl], axis=0), gradient[:,self.nr_parameters_ampl:].flatten()))

                H00 = hessian[:, :self.nr_parameters_ampl, :self.nr_parameters_ampl].sum(axis=0)
                H01 = np.concatenate([hessian[t, :self.nr_parameters_ampl, self.nr_parameters_ampl:] for t in range(self.nr_timeslots)], axis=1)
                H10 = np.concatenate([hessian[t, self.nr_parameters_ampl:, :self.nr_parameters_ampl] for t in range(self.nr_timeslots)], axis=0)
                H11 = scipy.linalg.block_diag(*[hessian[t, self.nr_parameters_ampl:, self.nr_parameters_ampl:] for t in range(self.nr_timeslots)])

                hessian = np.block([[H00, H01],[H10, H11]])
                hessian0 = hessian

                dx = np.dot(np.linalg.pinv(hessian, self.pinv_tol), gradient)
                norm_dx += np.linalg.norm(dx)**2

                if max_dx <  np.amax(abs(dx)) :
                    max_dx = np.amax(abs(dx))
                    i_max = i

                p0 = parameters[i].copy()

                parameters[i] += self.solver_update_gain*dx

                aterm_ampl = np.repeat(np.tensordot(parameters[i,:self.nr_parameters_ampl], self.Bampl, axes = ((0,), (0,)))[np.newaxis,:], self.nr_timeslots, axis=0)
                aterm_phase = np.exp(1j * np.tensordot(parameters[i,self.nr_parameters_ampl:].reshape((self.nr_timeslots, self.nr_parameters_phase)), self.Bphase, axes = ((1,), (0,))))

                aterms0 = aterms.copy()
                aterms[:, i] = aterm_ampl * aterm_phase

                timer1 += time.time()

            dresidual = previous_residual - residual_sum
            fractional_dresidual = dresidual / residual_sum

            print(max_dx, fractional_dresidual)


            previous_residual = residual_sum


            converged = (nr_iterations>1) and (fractional_dresidual < 1e-2)

            if converged:
                msg = "Converged after {nr_iterations} iterations - {max_dx}".format(nr_iterations=nr_iterations, max_dx=max_dx)
                print(msg)
                break

            if nr_iterations == 100:
                msg = "Did not converge after {nr_iterations} iterations - {max_dx}".format(nr_iterations=nr_iterations, max_dx=max_dx)
                print(msg)
                #figtitle.set_text(msg)
                break

        parameters_polynomial = parameters.copy()

        for i in range(self.nr_stations):
            parameters_polynomial[i,:self.nr_parameters_ampl] = np.dot(self.Tampl, parameters_polynomial[i,:self.nr_parameters_ampl])
            for j in range(self.nr_timeslots):
                parameters_polynomial[i,self.nr_parameters_ampl+j*self.nr_parameters_phase:self.nr_parameters_ampl+(j+1)*self.nr_parameters_phase] = \
                    np.dot(self.Tphase, parameters_polynomial[i,self.nr_parameters_ampl+j*self.nr_parameters_phase:self.nr_parameters_ampl+(j+1)*self.nr_parameters_phase])

        self.output_dataset.resize(self.output_dataset.shape[0]+1, axis=0)
        self.output_dataset[-1, :, :] = parameters_polynomial
        self.output_file.flush()


        timer += time.time()

        print("in {0} seconds".format(timer))
        print("in {0} seconds per iteration".format(timer/nr_iterations))
        print("spend {0} seconds in kernel".format(timer0))

    def update_info(self, dpinfo) :
        super().update_info(dpinfo)
        print("MyDPStep.update_info")
        self.info().set_need_vis_data()
        self.fetch_uvw = True
        self.fetch_weights = True
        self.frequencies = np.array(self.info().get_channel_frequencies()).astype(np.float32)
        self.nr_channels = self.frequencies.shape[0]
        self.antenna1 = np.array(self.info().get_antenna1())
        self.antenna2 = np.array(self.info().get_antenna2())
        self.nr_baselines = self.antenna1.shape[0]
        self.nr_stations = self.info().nantenna()
        self.time_start = self.info().startTime()
        self.time_step = self.info().timeInterval()

        self.baselines = np.zeros(shape=(self.nr_baselines), dtype=idg.baselinetype)
        self.baselines['station1'] = self.antenna1
        self.baselines['station2'] = self.antenna2



    def init_buffers(self):

        # Initialize empty buffers
        self.uvw          = np.zeros(shape=(self.nr_baselines, self.nr_timesteps),
                                dtype=idg.uvwtype)
        self.visibilities = np.zeros(shape=(self.nr_baselines, self.nr_timesteps, self.nr_channels,
                                    self.nr_correlations),
                                dtype=idg.visibilitiestype)
        self.weights      = np.zeros(shape=(self.nr_baselines, self.nr_timesteps, self.nr_channels,
                                    self.nr_correlations),
                                dtype=np.float32)

    def init_grid(self):
        h = fits.getheader(self.model_image)

        N0 = h["NAXIS1"]
        N = next_composite(int(N0 * self.padding))

        self.grid_size = N0
        self.padded_grid_size = N

        self.init_taper()

        self.img = np.zeros(shape=(self.nr_correlations, self.padded_grid_size, self.padded_grid_size), dtype=idg.gridtype)

        self.cell_size = abs(h["CDELT1"]) / 180 * np.pi
        self.image_size = self.padded_grid_size * self.cell_size
        d = fits.getdata(self.model_image)

        print("{0} non-zero pixels.".format(np.sum(d != 0.0)))

        self.img[0, (N-N0)//2:(N+N0)//2, (N-N0)//2:(N+N0)//2] = d[0,0,:,:]/np.outer(self.taper_grid0, self.taper_grid0)
        self.img[3, (N-N0)//2:(N+N0)//2, (N-N0)//2:(N+N0)//2] = d[0,0,:,:]/np.outer(self.taper_grid0, self.taper_grid0)

        self.grid = self.img.copy()
        self.proxy.transform(idg.ImageDomainToFourierDomain, self.grid)

    def init_taper(self):

        N0 = self.grid_size
        N = self.padded_grid_size
        # Initialize taper
        taper = idgwindow(self.subgrid_size, self.taper_support, self.padding)
        self.taper2 = np.outer(taper, taper).astype(np.float32)

        taper_ = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(taper)))
        taper_grid = np.zeros(N, dtype=np.complex128)
        taper_grid[(N-self.subgrid_size)//2:(N+self.subgrid_size)//2] = taper_ * np.exp(-1j*np.linspace(-np.pi/2, np.pi/2, self.subgrid_size, endpoint=False))
        taper_grid = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(taper_grid))).real*N/self.subgrid_size
        self.taper_grid0 = taper_grid[(N-N0)//2:(N+N0)//2]

    def init_idg_h5parm(self):


        self.output_filename = "parameters-0.01-{0}ampl-{1}phase.hdf5".format(self.nr_parameters_ampl, self.nr_parameters_phase)
        try:
            self.output_file = h5py.File(self.output_filename, "w")
            self.output_dataset = self.output_file.create_dataset("parameters", shape=(0,self.nr_stations, self.nr_parameters), maxshape=(None, self.nr_stations, self.nr_parameters), dtype=np.float32, chunks = True)
        except Exception as e:
            pass
            raise e


        self.output_dataset.attrs["nr_timeslots"]          = self.nr_timeslots
        self.output_dataset.attrs["nr_timesteps_per_slot"] = self.nr_timesteps_per_slot
        self.output_dataset.attrs["nr_parameters_ampl"]    = self.nr_parameters_ampl
        self.output_dataset.attrs["nr_parameters_phase"]   = self.nr_parameters_phase
        self.output_dataset.attrs["time_step"] = self.time_step
        self.output_dataset.attrs["time_start"] = self.time_start
        self.output_dataset.attrs["subgrid_size"] = self.subgrid_size
        self.output_dataset.attrs["image_size"] = self.image_size
        self.output_dataset.attrs["imagename"] = self.model_image
        self.output_dataset.attrs["solver_update_gain"] = self.solver_update_gain
        self.output_dataset.attrs["pinv_tol"] = self.pinv_tol

        self.output_file.flush()

    def init_basefunctions(self):
        B0 = np.ones((1, self.subgrid_size, self.subgrid_size,1))

        x = np.linspace(-0.5, 0.5, self.subgrid_size)

        B1,B2 = np.meshgrid(x,x)

        B1 = B1[np.newaxis, :, :, np.newaxis]
        B2 = B2[np.newaxis, :, :, np.newaxis]
        B3 = B1*B1
        B4 = B2*B2
        B5 = B1*B2
        B6 = B1*B1*B1
        B7 = B1*B1*B2
        B8 = B1*B2*B2
        B9 = B2*B2*B2
        B10 = B1*B1*B1*B1
        B11= B1*B1*B1*B2
        B12 = B1*B1*B2*B2
        B13 = B1*B2*B2*B2
        B14 = B2*B2*B2*B2
        B15 = B1*B1*B1*B1*B1
        B16 = B1*B1*B1*B1*B2
        B17 = B1*B1*B1*B2*B2
        B18 = B1*B1*B2*B2*B2
        B19 = B1*B2*B2*B2*B2
        B20 = B2*B2*B2*B2*B2

        base_polynomial = np.concatenate((B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20))

        base_polynomial_ampl = base_polynomial[:self.nr_parameters_ampl].reshape((-1, self.subgrid_size*self.subgrid_size)).T
        U,S,V, = np.linalg.svd(base_polynomial_ampl)
        base_orthonormal_ampl = U[:, :self.nr_parameters_ampl]
        self.Tampl = np.dot(np.linalg.pinv(base_polynomial_ampl), base_orthonormal_ampl)
        base_orthonormal_ampl = base_orthonormal_ampl.T.reshape((-1, self.subgrid_size,self.subgrid_size, 1))

        base_polynomial_phase = base_polynomial[:self.nr_parameters_phase].reshape((-1, self.subgrid_size*self.subgrid_size)).T
        U,S,V, = np.linalg.svd(base_polynomial_phase)
        base_orthonormal_phase = U[:, :self.nr_parameters_phase]
        self.Tphase = np.dot(np.linalg.pinv(base_polynomial_phase), base_orthonormal_phase)
        base_orthonormal_phase = base_orthonormal_phase.T.reshape((-1, self.subgrid_size,self.subgrid_size, 1))

        self.Bampl = np.kron(base_orthonormal_ampl, np.array([1.0, 0.0, 0.0, 1.0]))
        self.Bphase = np.kron(base_orthonormal_phase, np.array([1.0, 0.0, 0.0, 1.0]))

import pydppp
import numpy as np

class MyDPStep(pydppp.DPStep):

    def __init__(self,parset, prefix):
        super().__init__()
        self.parset = parset
        self.dpbuffers = []
        print(dir(parset))
        print("prefix = ", prefix)
        print(parset.getString(prefix + "type"))

    def show(self) :
        print("MyDPStep.show")

    def process(self, dpbuffer) :
        pass
        print("process", np.array(dpbuffer.get_flags(), copy=False).shape)
        #print("process", self.get_count(), np.array(self.get_uvw(), copy=False).shape)
        self.dpbuffers.append(dpbuffer)
        if len(self.dpbuffers) == 10:
            self.calibrate()
            for dpbuffer in self.dpbuffers:
                pass
                self.process_next_step(dpbuffer)

            self.dpbuffers = []


    def calibrate(self):
        pass

    def update_info(self, dpinfo) :
        super().update_info(dpinfo)
        print("MyDPStep.update_info")
        self.info().set_need_vis_data()
        self.fetch_uvw = True
        self.fetch_weigths = True

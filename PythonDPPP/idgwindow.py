import numpy

def idgwindow(N, W, padding, offset = 0.5, l_range = None):
  """

  """
  
  l_range_inner = numpy.linspace(-(1/padding)/2,(1/padding)/2, N*16+1)

  vl = (numpy.arange(N)-N/2+offset)/N
  vu = numpy.arange(N) - N/2 + offset
  Q = numpy.sinc((N-W+1)*(vl[:,numpy.newaxis] - vl[numpy.newaxis,:]))
  
  B = []
  RR = []
  for l in l_range_inner:
      d = numpy.mean(numpy.exp(2*numpy.pi*1j*vu[numpy.newaxis,:] * (vl[:, numpy.newaxis]-l)), axis=1).real
      D = d[:,numpy.newaxis] * d[numpy.newaxis,:]
      b_avg = numpy.sinc((N-W+1)*(l-vl))
      B.append(b_avg*d)
      S = b_avg[:,numpy.newaxis] * b_avg[numpy.newaxis,:]
      RR.append(D*(Q - S))
  B = numpy.array(B)
  RR = numpy.array(RR)
      
  taper = numpy.ones(len(l_range_inner))

  for q in range(10):
      R = numpy.sum((RR*1/taper[:,numpy.newaxis, numpy.newaxis]**2), axis = 0)
      R1 = R[:,:(N//2)] + R[:,:(N//2)-1:-1]
      R2 = R1[:(N//2),:] + R1[:(N//2)-1:-1,:]
      U, S1, V = numpy.linalg.svd(R2)
      a = numpy.abs(numpy.concatenate([U[:,-1],U[::-1,-1]]))
      taper = numpy.dot(B,a)

  if l_range is None:
    return a
  else:
    B = []
    RR = []
    for l in l_range:
        d = numpy.mean(numpy.exp(2*numpy.pi*1j*vu[numpy.newaxis,:] * (vl[:, numpy.newaxis]-l)), axis=1).real
        D = d[:,numpy.newaxis] * d[numpy.newaxis,:]
        b_avg = numpy.sinc((N-W+1)*(l-vl))
        B.append(b_avg*d)
        S = b_avg[:,numpy.newaxis] * b_avg[numpy.newaxis,:]
        RR.append(D*(Q - S))
    B = numpy.array(B)
    RR = numpy.array(RR)
    
    return a, B, RR


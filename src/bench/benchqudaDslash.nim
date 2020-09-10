#[
import base
import layout
import gauge
import strUtils
import gauge/fat7l
import quda/quda_milc_interface
import quda/quda
import quda/enum_quda
import quda/qudaWrapper
import physics/qcdTypes
import gauge/staples
import maths

export PerfInfo

type
  D4ColorMatrix = array[4, DColorMatrix]
  D4LatticeColorMatrix = Field[1, D4ColorMatrix]


proc benchmark( g_in: any, r_in,t_in:Field, info: var PerfInfo) =
  let lo = r_in.l.qudaSetup
  var
    gauge_param: ptr QudaGaugeParam
    g : D4LatticeColorMatrix
    link : pointer = g.s.data
    t, r: DLatticeColorVector
  t.new lo
  r.new lo
  g.new lo
  var
    invargs: QudaInvertArgs_t
    precision = 2
    srcGpu: pointer = t.s.data
    destGpu: pointer = r.s.data
    fatlink: pointer = g.s.data
    longlink: pointer = nil
  invargs.evenodd = QUDA_EVEN_PARITY


  threads:

    for i in r.sites:
      var cv: array[4,cint]
      r_in.l.coord(cv,(r_in.l.myRank,i))
      let ri = lo.rankIndex(cv)
      forO a, 0, 2:
        t[ri.index][a].re := t_in{i}[a].re
        t[ri.index][a].im := t_in{i}[a].im
    for i in r.sites:
      var cv: array[4,cint]
      r_in.l.coord(cv,(r_in.l.myRank,i))
      let ri = lo.rankIndex(cv)
      forO a, 0, 2:
        r[ri.index][a].re := r_in{i}[a].re
        r[ri.index][a].im := r_in{i}[a].im

    if g_in.len == 4: # plain staggered
      for i in t.l.sites:
        var cv: array[4,cint]
        r_in.l.coord(cv,(r_in.l.myRank,i))
        let ri = lo.rankIndex(cv)
        forO mu, 0, 3:
          forO a, 0, 2:
            forO b, 0, 2:
              g[ri.index][mu][a,b].re := g_in[mu]{i}[a,b].re
              g[ri.index][mu][a,b].im := g_in[mu]{i}[a,b].im
    else:
      echo "unsupported g_in.len: ", g_in.len
      quit(-1)

    var dummy : cint = -1
    qudaDslash(precision.cint, precision.cint,invargs,fatlink,longlink,srcGpu,destGpu,addr(dummy))

    threads:
      for i in r.sites:
        var cv: array[4,cint]
        r_in.l.coord(cv,(r_in.l.myRank,i))
        let ri = lo.rankIndex(cv)
        forO a, 0, 2:
          r_in{i}[a].re := r[ri.index][a].re
          r_in{i}[a].im := r[ri.index][a].im


when isMainModule:
  import qex
  import physics/qcdTypes
  import gauge
  qexInit()
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  #for mu in 0..<g.len: g[mu] := 1
  var
    src = lo.ColorVector()
    dest = lo.ColorVector()
    rng = RngMilc6.newRNGField lo
    info : PerfInfo
  g.random rng
  threads:
    g.setBC
    g.stagPhase
    dest := 0
    src.z4 rng

  benchmark(g,dest,src,info)


  qexFinalize()
]#

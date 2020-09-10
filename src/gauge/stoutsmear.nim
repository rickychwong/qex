import base
import layout
import gauge
import strUtils
import fat7l
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


type StoutCoefs* = object
  n*: cuint
  rho*: cdouble
  layout* : Layout[4]

proc `$`*(c: StoutCoefs): string =
  result = "Stout{\n"
  result &= "  n: " & $c.n & "\n"
  result &= "  rho: " & $c.rho & "\n"
  result &= "}"

proc smearGetForce*[G](coef: StoutCoefs, gf: G, fl:G,
            info: var PerfInfo):auto =
  ## Note that the resulting proc, smearedForce, holds a reference to the input gauge gf.
  ## The correctness of the algorithm depends on gf remaining the same.
  ## On the contrary, any changes to the smeared gauge fl would have no effects to the force calculation.
  
  proc smearing(result, g: G)=
    var
      cs = startCornerShifts(g)
      (stf,stu,ss) = makeStaples(g,cs)
    let nd = coef.layout.nDim
    threads:
      #for ir in coef.layout.sites:
      for ir in g[0]:
        var tmp : type(load1(g[0][0]))
        tmp := 0
        for mu in 0..<nd:
          for nu in 0..<nd:
            if mu!=nu:
              if isLocal(ss[mu][nu],ir):
                var tmp2 : type(load1(g[0][0])) 
                tmp2 := 0
                var tmp3 : type(load1(g[0][0])) 
                tmp3 := 0
                mul(tmp2,g[mu][ir],stf[mu,nu][ir].adj)
                add(tmp,tmp,tmp2)
                
                localSB(ss[mu][nu], ir, assign(tmp3,it), stu[mu,nu][ix])
                mul(tmp2,g[mu][ir],tmp3.adj)
                add(tmp,tmp,tmp2)
          for nu in 0..<nd:
            if mu!=nu:
              var needBoundary = false
              boundaryWaitSB(ss[mu][nu]): needBoundary = true
              if needBoundary:
                for ir in coef.layout:
                  if not isLocal(ss[mu][nu],ir):
                    var tmp2 : type(load1(g[0][0])) 
                    tmp2 := 0
                    var tmp3 : type(load1(g[0][0])) 
                    tmp3 := 0
                    getSB(ss[mu][nu], ir, assign(tmp3,it), stu[mu,nu][ix])
                    mul(tmp2,g[mu][ir],tmp3.adj)
                    add(tmp,tmp,tmp2)


          #result[mu][i] := exp( -coef.rho * tmp ) * g[mu][i]
          mul( result[mu][ir] , exp( -coef.rho * tmp ) , g[mu][ir] )

  smearing(fl, gf)
  for i in 0..<coef.n:
    smearing(fl,fl)

  #[ 
    tic()
    let lo = coef.layout.qudaSetup
    var
      gauge_param: ptr QudaGaugeParam
      g : D4LatticeColorMatrix
      link : pointer = g.s.data

    g.new lo
    threads:
      if gf.len == 4: # plain staggered
        link = nil
        for i in coef.layout.sites:
            var cv: array[4,cint]
            lo.coord(cv,(lo.myRank,i))
            let ri = lo.rankIndex(cv)
            forO mu, 0, 3:
              forO a, 0, 2:
                forO b, 0, 2:
                  g[ri.index][mu][a,b].re := gf[mu]{i}[a,b].re
                  g[ri.index][mu][a,b].im := gf[mu]{i}[a,b].im
      else:
          echo "unsupported gf.len: ", gf.len
          quit(-1)

    loadGaugeQuda(link, gauge_param)
    performSTOUTnStep(coef.n, coef.rho )
    saveGaugeQuda(link, gauge_param)

    threads:
      if fl.len == 4: # plain staggered
        link = nil
        for i in coef.layout.sites:
            var cv: array[4,cint]
            lo.coord(cv,(lo.myRank,i))
            let ri = lo.rankIndex(cv)
            forO mu, 0, 3:
              forO a, 0, 2:
                forO b, 0, 2:
                  fl[mu]{i}[a,b].re := g[ri.index][mu][a,b].re 
                  fl[mu]{i}[a,b].im := g[ri.index][mu][a,b].im
      else:
          echo "unsupported fl.len: ", fl.len
          quit(-1)
  
  ]#

  proc smearedForce(f,chain:G) = 
    tic("smearedF")
  smearedForce




proc smearPriv[G](coef: StoutCoefs, gf: G, fl: G, info: var PerfInfo) {.codegenDecl:
    "__attribute__((noinline)) $# $#$#".} =
  # Avoid inlining or other compiler optimizations
  # in order to guarantee the change of the stack pointer,
  # such that Nim's GC is able to collect the memory.
  {.emit: "asm (\"\");".}
  var f = coef.smearGetForce(gf, fl, info)
  f = nil

proc smear*[G](coef: StoutCoefs, gf: G, fl: G, info: var PerfInfo) =
  ## Try our best to release memory here.
  ## Sometimes it still requires a GC after this function returns.
  coef.smearPriv(gf, fl, info)
  qexGC()

proc smear*(c: StoutCoefs, g: any, fl: any) =
  var info: PerfInfo
  c.smear(g, fl, info)

when isMainModule:
  import qex
  import physics/qcdTypes
  import gauge
  qexInit()
  #var defaultGaugeFile = "l88.scidac"
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  #for mu in 0..<g.len: g[mu] := 1
  #g.random

  var info: PerfInfo
  var coef: StoutCoefs
  coef.n = 2
  coef.rho = 0.15

  echo coef

  var fl = lo.newGauge()

  template disp(g: typed) =
    let p = g.plaq
    let sp = 2.0*(p[0]+p[1]+p[2])
    let tp = 2.0*(p[3]+p[4]+p[5])
    echo p
    echo sp
    echo tp
    echo trace(g[0])

  disp g
  coef.smear(g, fl, info)
  disp fl
  qexFinalize()

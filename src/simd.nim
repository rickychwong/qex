#import simdGcc
#export simdGcc

when defined(SSE) or defined(AVX) or defined(AVX512):
  import simd/simdX86
  export simdX86
elif defined(QPX):
  import simd/simdQpx
  export simdQpx

#import simd/simdGeneric
#export simdGeneric


proc toSingle*(x: SimdD4): SimdS4 {.inline,noInit.} =
  for i in 0..<4:
    result[i] = x[i]
proc toSingle*(x: SimdD8): SimdS8 {.inline,noInit.} =
  for i in 0..<8:
    result[i] = x[i]

when declared(SimdS4):
  proc toDouble*(x: SimdS4): SimdD4 {.inline,noInit.} =
    for i in 0..<4:
      result[i] = x[i]
  proc inorm2*(r:var SimdD4; x:SimdS4) {.inline.} =
    let y = toDouble(x)
    inorm2(r, y)

template assign*(r: array[4,float32], x: SimdD4): untyped =
  assign(r, toSingle(x))
template assign*(r: SimdS4, x: SimdD4): untyped =
  r = toSingle(x)
template assign*(r: SimdS8, x: SimdD8): untyped =
  r = toSingle(x)
template assign*(r: SimdD4, x: SimdS4): untyped =
  r = toDouble(x)
template assign*(r: SimdD8, x: SimdS8): untyped =
  r = toDouble(x)

template isub*(r: SimdD8, x: SimdS8): untyped =
  isub(r, toDouble(x))
template imadd*(r: SimdD8, x: SimdD8, y: SimdS8): untyped =
  imadd(r, x, toDouble(y))
template imsub*(r: SimdD8, x: SimdD8, y: SimdS8): untyped =
  imsub(r, x, toDouble(y))

template eval*(x: SimdD8): auto = x

when declared(SimdS4):
  converter promote*(x: SimdS4): SimdD4 {.inline,noInit.} =
    assign(result, x)
  template toSingleImpl*(x: SimdS4): untyped = x
  template toSingleImpl*(x: SimdD4): untyped = toSingle(x)
  template toDoubleImpl*(x: SimdS4): untyped = toDouble(x)
  template toDoubleImpl*(x: SimdD4): untyped = x

when declared(SimdS8):
  converter promote*(x: SimdS8): SimdD8 =
    assign(result, x)
  template toSingleImpl*(x: SimdS8): untyped = x
  template toSingleImpl*(x: SimdD8): untyped = toSingle(x)
  template toDoubleImpl*(x: SimdS8): untyped = toDouble(x)
  template toDoubleImpl*(x: SimdD8): untyped = x

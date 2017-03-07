import osPaths
import strUtils

var params = newSeq[string](0)
template set(key,val: string) =
  params.insert val
  params.insert key
template `~`(key,val: untyped) =
  set(astToStr(key), val)

let nim = paramStr(0)
NIM ~ nim

let script = paramStr(1)
#let qexdir = thisDir()
let (qexdir, _, _) = splitFile(script)
echo "Using QEXDIR ", qexdir
set "QEXDIR", qexdir

let home = getHomeDir()

var qmpdir = home / "lqcd/install/qmp"
echo "Using QMPDIR ", qmpdir
set "QMPDIR", qmpdir

var qiodir = home / "lqcd/install/qio"
echo "Using QIODIR ", qiodir
set "QIODIR", qiodir

var machine = ""
CC ~ "mpicc"
LD ~ "$(CC)"
CC_TYPE ~ "gcc"
CFLAGS_ALWAYS ~ "-Wall -std=gnu99 -march=native -ldl"
CFLAGS_DEBUG ~ "-g3 -O0"
CFLAGS_SPEED ~ "-g -O3"
VERBOSITY ~ "1"
SIMD ~ ""
VLEN ~ "1"

if dirExists "/bgsys":
  machine = "Blue Gene/Q"
  CC ~ qexdir / "mpixlc2"
  CFLAGS_ALWAYS ~ ""
  CFLAGS_DEBUG ~ "-g3 -O0"
  CFLAGS_SPEED ~ "-O3"
  VERBOSITY ~ "3"
  SIMD ~ "QPX"

let uname = staticExec "uname"

if machine=="" and uname=="Darwin":
  machine = "macOS"
  SIMD ~ "SSE,AVX"
  VLEN ~ "8"

if machine=="" and fileExists "/proc/cpuinfo":
  # assume compute nodes are same as build nodes
  machine = "linux"
  SIMD ~ "SSE,AVX"
  VLEN ~ "8"

# cray/modules
# KNL
# linux (/proc/cpuinfo)
# check on linux/mac/vesta/cooley/theta
# gcc/icc opt level

var f = readFile(qexdir / "Makefile.in")
f = replace(f, "$", "!DOLLAR!")
f = replace(f, "#", "!HASH!")
f = replace(f, "@@", "$")
f = f % params
f = replace(f, "!HASH!", "#")
f = replace(f, "!DOLLAR!", "$")

#echo f
writeFile("Makefile", f)

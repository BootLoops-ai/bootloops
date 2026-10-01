// ginac_gpl.cpp — tiny CLI: G(a₁,…,aₙ; z) at arbitrary precision via GiNaC.
//
// License note: this source file is MIT (part of BootLoops), but it links the
// GiNaC and CLN libraries (GPL-2+), so any compiled binary is a GPL combined
// work.  The binary is built locally on first use by ginac_bridge.jl and is
// not distributed; GPLEval.jl only talks to it as a subprocess.
//
//   ginac_gpl <digits> <re_z> <im_z> <re_a1> <im_a1> [<re_a2> <im_a2> ...]
//
// Prints "Re\nIm\n" of G(a;z) at <digits> decimal digits (CLN). Letters and z
// are parsed as exact rationals if they look like p/q, otherwise as decimal
// strings (CLN reads them at the requested precision).
//
// Compile:  g++ -O2 -std=c++17 -o ginac_gpl ginac_gpl.cpp $(pkg-config --cflags --libs ginac)

#include <ginac/ginac.h>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace GiNaC;

static numeric parse_num(const char* s) {
    return numeric(s);   // CLN parses "p/q" exactly, decimals at current Digits
}

int main(int argc, char** argv) {
    if (argc < 4 || ((argc - 4) % 2) != 0) {
        std::cerr << "usage: ginac_gpl digits re_z im_z [re_a im_a]...\n";
        return 1;
    }
    Digits = std::atol(argv[1]) + 10;
    numeric zr = parse_num(argv[2]), zi = parse_num(argv[3]);
    ex z = zr + I * zi;
    lst a;
    for (int i = 4; i + 1 < argc; i += 2) {
        numeric ar = parse_num(argv[i]), ai = parse_num(argv[i + 1]);
        a.append(ar + I * ai);
    }
    ex g = (a.nops() == 0) ? ex(1) : G(a, z);
    ex v = g.evalf();
    if (!is_a<numeric>(v)) { std::cerr << "non-numeric result: " << v << "\n"; return 2; }
    numeric n = ex_to<numeric>(v);
    std::cout.precision(Digits);
    std::cout << n.real() << "\n" << n.imag() << "\n";
    return 0;
}

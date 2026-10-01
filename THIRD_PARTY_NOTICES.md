# Third-party notices

DADA remains GPL-3.0-or-later; the canonical LICENSE is unchanged.

## CoolProp 8.0.0 — MIT

`src/dada_solver/numerical_primitives.py` adapts the dilute air coefficients
from `dev/fluids/Air.json` and dilute air/helium formulas from
`src/Backends/Helmholtz/TransportRoutines.cpp` at tag `v8.0.0`:
https://github.com/CoolProp/CoolProp/tree/v8.0.0

Only these formulas/coefficients are incorporated. CoolProp is an optional
reference-generation dependency, not a runtime dependency. No third-party
components of CoolProp are vendored.

MIT License

Copyright (c) 2012-2018 Ian H. Bell and other CoolProp developers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

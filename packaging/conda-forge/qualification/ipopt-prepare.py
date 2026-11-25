"""Only increase reporting precision and explicitly select MUMPS in archive HS071."""
from pathlib import Path
p = Path('ipopt-source/examples/hs071_cpp/hs071_main.cpp')
s = p.read_text()
if 'std::setprecision(17)' in s:
    raise SystemExit(0)
s = s.replace('#include "IpIpoptApplication.hpp"', '#include "IpIpoptApplication.hpp"\n#include <iomanip>')
s = s.replace('  // Create a new instance of your nlp', '  std::cout << std::setprecision(17);\n  // Create a new instance of your nlp', 1)
s = s.replace('   app->Options()->SetNumericValue("tol", 3.82e-6);', '   app->Options()->SetNumericValue("tol", 1e-7);\n  app->Options()->SetStringValue("linear_solver", "mumps");\n  app->Options()->SetStringValue("option_file_name", "");')
p.write_text(s)

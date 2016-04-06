#!/usr/bin/env python
import os
import sys
import imp


tool = sys.argv[1]
methgo_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = 'scripts/{0}/{0}.py'.format(tool)
filename = os.path.join(methgo_dir, module_dir)
module = imp.load_source(tool, filename)
del sys.argv[0]
module.main()

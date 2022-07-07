from sympy import *
init_printing(use_latex=True)
from pprint import pp

x, y = symbols('x y', real=True)


def complsq(expr, var):
  
  # expr is of the form ax^2 + bx + c, and var == x
  
  error_msg = '''
  Expression must be of the form ax^2 + bx + c.
  The input was ''' + latex(expr)
  
  # var may be any symbol, but let's call it x here
  x = var
  
  # Used to match the coefficients
  a, b, c = [ \
    Wild(w, exclude=[x, y]) for w in 'abc' \
  ]
  
  # Match
  m = expr.match(a*x**2 + b*x + c)
  
  # No match: bad expression
  if m is None:
    raise ValueError(error_msg)
  
  # Get coefficients
  A, B, C = [w.xreplace(m) for w in [a, b, c]]
  
  # Expression must be of degree 2
  if A == 0:
    raise ValueError(error_msg)
  
  # Will return A (x + B/2A)^2 + (4AC - B^2)/4A
  second = sympify(B) / (2*A)
  third = sympify(4*A*C - B**2) / (4*A)
  
  return A * (var + second)**2 + third


def classify(eq):
  
  error_msg = '''
  Equation must be of the form Ax^2 + Cy^2 + Dx + Ey + F.
  The input was ''' + latex(eq)
  
  # Used to match the coefficients of the equation
  a, b, c, d, e, f = [ \
    Wild(w, exclude=[x, y]) for w in 'abcdef' \
  ]
  
  # Match
  m = eq.match(a*x**2 + b*x*y + c*y**2 + d*x + e*y + f)
  
  # No match: bad equation
  if m is None:
    raise ValueError(error_msg)
  
  # Get coefficients
  A, B, C, D, E, F = [w.xreplace(m) for w in [a, b, c, d, e, f]]
  
  # Must be degree 2 in at least one of the vars
  # No rotated conics for now
  if B != 0 or (A == 0 and C == 0):
    raise ValueError(error_msg)
  
  # Collect everything in a dict
  retval = { \
    'eq_general': Eq(eq, 0), \
    'A': A, \
    'B': B, \
    'C': C, \
    'D': D, \
    'E': E, \
    'F': F \
  }
  
  # Discriminants
  # Non-homogeneous form
  M = Matrix([[A, B/2], [B/2, C]])
  detM = M.det()
  
  # Homogeneous form
  Q = Matrix([[A, B/2, D/2], [B/2, C, E/2], [D/2, E/2, F]])
  detQ = Q.det()
  
  # Degenerate?
  is_degenerate = detQ == 0
  
  # We don't know which form yet
  degenerate_form = None
  
  # Type
  if detM > 0:
    # Ellipse or circle
    # Must check if empty set
    temp = A * (E**2 / (4*C) + D**2 / (4*A) - F)
    if temp < 0:
      conic_type = 'ellipse'
      is_degenerate = True
      degenerate_form = 'empty set'
    else:
      if A == C:
        conic_type = 'circle'
      else:
        conic_type = 'ellipse'
  elif detM == 0:
    conic_type = 'parabola'
  else:
    conic_type = 'hyperbola'
    
  # Degenerate form (but not empty-set ellipse, which has already been
  # identified above)
  if is_degenerate and degenerate_form is None:
    
    if conic_type == 'hyperbola':
      degenerate_form = 'intersecting lines'
      
    elif conic_type == 'parabola':
      temp = D**2 + E**2 - 4 * F * (A + C)
      if temp < 0:
        degenerate_form = 'empty set'
      elif temp == 0:
        degenerate_form = 'coincident lines'
      else: # temp > 0
        degenerate_form = 'parallel lines'
        
    elif detM > 0:
      conic_type = 'ellipse'
      degenerate_form = 'point'
      
  # Save info in dict
  retval['conic_type'] = conic_type
  retval['is_degenerate'] = is_degenerate
  if is_degenerate:
    retval['degenerate_form'] = degenerate_form
  
  return retval


def find_canonical_eq(eq):
  
  conic = classify(eq)
  conic_type = conic['conic_type']
  
  if conic_type == 'parabola':
    do_parabola(conic)
  elif conic_type in ['ellipse', 'circle']:
    do_ellipse(conic)
  elif conic_type == 'hyperbola':
    do_hyperbola(conic)
  else:
    raise ValueError('Impossible conic type!')
  
  return conic


def do_parabola(conic):
  
  A, C, D, E, F = [conic[c] for c in 'ACDEF']
  
  # x part
  eqx = A*x**2 + D*x
  
  # y part
  eqy = C*y**2 + E*y
  
  # Used to build the new equation
  a, b, c = [Wild(w, exclude=[x, y]) for w in 'abc']
  
  if C == 0: # This is a parabola with x^2
    
    # Build new equation
    if not conic['is_degenerate']:  # Nondegenerate case
      
      # Complete the square only for the x part
      eqx = complsq(eqx, x)
      
      # Get constants in the completed square
      m = eqx.match(a * (x + b)**2 + c)
          
      neweq = Eq(y + (F + c)/E, -a * (x + b)**2 / E).xreplace(m)
      conic['eq_canonical'] = [neweq]
      
    else: # Degenerate cases (except empty set)
      
      if conic['degenerate_form'] == 'coincident lines':
        conic['eq_canonical'] = [Eq(x, -D/(2*A))]
      elif conic['degenerate_form'] == 'parallel lines':
        conic['eq_canonical'] = [ \
          Eq(x, (-D + sqrt(D**2 - 4*A*F))/(2*A)), \
          Eq(x, (-D - sqrt(D**2 - 4*A*F))/(2*A)) \
        ]
      
  else: # This is a parabola with y^2
    
    # Build new equation
    if not conic['is_degenerate']: # Nondegenerate case
          
      # Complete the square only for the y part
      eqy = complsq(eqy, y)
      
      # Get constants in the completed square
      m = eqy.match(a * (y + b)**2 + c)
      
      neweq = Eq(x + (F + c)/D, -a * (y + b)**2 / D).xreplace(m)
      conic['eq_canonical'] = [neweq]
      
    else: # Degenerate cases (except empty set)
      
      if conic['degenerate_form'] == 'coincident lines':
        conic['eq_canonical'] = [Eq(y, -E/(2*C))]
      elif conic['degenerate_form'] == 'parallel lines':
        conic['eq_canonical'] = [ \
          Eq(y, (-E + sqrt(E**2 - 4*C*F))/(2*C)), \
          Eq(y, (-E - sqrt(E**2 - 4*C*F))/(2*C)) \
        ]


def do_ellipse(conic):
  
  A, C, D, E, F = [conic[c] for c in 'ACDEF']
  
  # x part
  eqx = A*x**2 + D*x
  
  # y part
  eqy = C*y**2 + E*y
  
  # Used to build the new equation 
  a, b, r = [Wild(w, exclude=[x, y]) for w in 'abr']
  qx = Wild('qx', exclude=[y])
  qy = Wild('qy', exclude=[x])
  
  # Build a(x - ...)^2 + b(y - ...)^2 = r
  neweq = complsq(eqx, x) + complsq(eqy, y) + F
  m = neweq.match(a * qx**2 + b * qy**2 + r)
  
  aa = a.xreplace(m)
  bb = b.xreplace(m)
  rr = r.xreplace(m)
  qqx = qx.xreplace(m)
  qqy = qy.xreplace(m)
  
  # Make sure lhs will be positive
  if aa < 0 and bb < 0:
    aa, bb, rr = -aa, -bb, -rr
  
  neweq = Eq( \
    aa * qqx**2 + bb * qqy**2, -rr, \
    # Must keep unevaluated, or empty sets will simplify to False:
    evaluate=False \
  )
  
  # If nondegenerate ellipse, make rhs 1
  # And push constants to denominators
  if conic['conic_type'] == 'ellipse' and -rr > 0:
    
    neweq = Eq(neweq.lhs / -rr, 1)
    
    mm = neweq.lhs.match(qx**2 / a + qy**2 / b)
    aa = a.xreplace(mm)
    bb = b.xreplace(mm)
    qqx = qx.xreplace(mm)
    qqy = qy.xreplace(mm)
    
    newlhs = sympify( \
      'qqx**2 / aa + qqy**2 / bb', \
      locals={'qqx': qqx, 'aa': aa, 'qqy': qqy, 'bb': bb}, \
      evaluate=False \
    )
    
    neweq = Eq(newlhs, 1)
    
  elif conic['conic_type'] == 'circle':
    # If circle, make aa = bb = 1
    neweq = Eq(neweq.lhs / aa, -rr / aa)
    
  # Store in dict
  conic['eq_canonical'] = [neweq]
  
  # If point, find coordinates and append to field in dict
  if conic['is_degenerate'] and conic['degenerate_form'] == 'point':
    p = Point(solve(qqx)[0], solve(qqy)[0])
    conic['eq_canonical'].append(p)


def do_hyperbola(conic):
  
  A, C, D, E, F = [conic[c] for c in 'ACDEF']
  
  # x part
  eqx = A*x**2 + D*x
  
  # y part
  eqy = C*y**2 + E*y
  
  # Used to build the new equation 
  a, b, r = [Wild(w, exclude=[x, y]) for w in 'abr']
  qx = Wild('qx', exclude=[y])
  qy = Wild('qy', exclude=[x])
  
  # Build a(x - ...)^2 + b(y - ...)^2 = r
  neweq = complsq(eqx, x) + complsq(eqy, y) + F
  m = neweq.match(a * qx**2 + b * qy**2 + r)
  
  aa = a.xreplace(m)
  bb = b.xreplace(m)
  rr = r.xreplace(m)
  qqx = qx.xreplace(m)
  qqy = qy.xreplace(m)
  
  # Make sure rhs will be positive
  if rr > 0:
    aa, bb, rr = -aa, -bb, -rr
  
  neweq = Eq(aa * qqx**2 + bb * qqy**2, -rr)
  
  # If nondegenerate hyperbola, make rhs 1
  # And push constants to denominators
  if not conic['is_degenerate']:
    
    neweq = Eq(neweq.lhs / -rr, 1)
    
    mm = neweq.lhs.match(qx**2 / a + qy**2 / b)
    aa = a.xreplace(mm)
    bb = b.xreplace(mm)
    qqx = qx.xreplace(mm)
    qqy = qy.xreplace(mm)
    
    newlhs = sympify( \
      'qqx**2 / aa + qqy**2 / bb', \
      locals={'qqx': qqx, 'aa': aa, 'qqy': qqy, 'bb': bb}, \
      evaluate=False \
    )
    
    neweq = Eq(newlhs, 1)
    
    # Store in dict
    conic['eq_canonical'] = [neweq]
    
  else:
    
    # If degenerate, compute intersecting lines
    if aa < 0:
      l1 = sqrt(bb) * qqy - sqrt(-aa) * qqx
      l2 = sqrt(bb) * qqy + sqrt(-aa) * qqx
    else:
      l1 = sqrt(aa) * qqx - sqrt(-bb) * qqy
      l2 = sqrt(aa) * qqx + sqrt(-bb) * qqy
    
    conic['eq_canonical'] = [ \
      simplify(Eq(l1, 0)), simplify(Eq(l2, 0)) \
    ]

%Calculates the integration of the quartic interpolation function
% x^4 - 2 x^2 + 1, over the interval a to b.

full = 0

if full
  syms x_1 p_1 m_1 x_2 p_2 m_2 x_3 p_3 m_3 t

  f = ((x_1 - (p_1 + m_1 * t)) ^ 4 - 2 * (x_1 - (p_1 + m_1 * t)) ^ 2 + 1) * ...
      ((x_2 - (p_2 + m_2 * t)) ^ 4 - 2 * (x_2 - (p_2 + m_2 * t)) ^ 2 + 1) * ...
      ((x_3 - (p_3 + m_3 * t)) ^ 4 - 2 * (x_3 - (p_3 + m_3 * t)) ^ 2 + 1);

  f_part = ((x_1 - (p_1 + m_1 * t)) ^ 4 - 2 * (x_1 - (p_1 + m_1 * t)) ^ 2 + 1);  

  disp('f = ');

  pretty(f);

  f_int = int(f, t);

  % pretty(f_int);

  % disp('f_int_simplify = ');
  % 
  f_int = simplify(f_int);
  % 
  % pretty(f_int);
else
  
  syms a_1 b_1 c_1 d_1 e_1 a_2 b_2 c_2 d_2 e_2 a_3 b_3 c_3 d_3 e_3 t
  
  f = (a_1 * t^4 + b_1 * t^3 + c_1 * t^2 + d_1 * t + e_1) *...
      (a_2 * t^4 + b_2 * t^3 + c_2 * t^2 + d_2 * t + e_2) *...
      (a_3 * t^4 + b_3 * t^3 + c_3 * t^2 + d_3 * t + e_3);
    
  pretty(f);
  
  f = collect(expand(f),t)
  
  
end

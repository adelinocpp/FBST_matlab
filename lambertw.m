function w = lambertw_aps(b,z)
%LAMBERTW  Lambert W-Function.
%   W = LAMBERTW(Z) computes the principal value of the Lambert 
%   W-Function, the solution of Z = W*exp(W).  Z may be a 
%   complex scalar or array.  For real Z, the result is real on
%   the principal branch for Z >= -1/e.
%
%   W = LAMBERTW(B,Z) specifies which branch of the Lambert 
%   W-Function to compute.  If Z is an array, B may be either an
%   integer array of the same size as Z or an integer scalar.
%   If Z is a scalar, B may be an array of any size.
%
%   The algorithm uses series approximations as initializations
%   and Halley's method as developed in Corless, Gonnet, Hare,
%   Jeffrey, Knuth, "On the Lambert W Function", Advances in
%   Computational Mathematics, volume 5, 1996, pp. 329-359.

% Pascal Getreuer 2005-2006
% Modified by Didier Clamond, 2005

if nargin == 1
   z = b;
   b = 0;
end

% Use asymptotic expansion w = log(z) - log(log(z)) for most z
tmp = log(z + (z == 0)) + i*b*6.28318530717958648;
w = tmp - log(tmp + (tmp == 0));

% For b = 0, use a series expansion when close to the branch point
k = find(b == 0 & abs(z + 0.3678794411714423216) <= 1.5);
tmp = sqrt(5.43656365691809047*z + 2) - 1 + i*b*6.28318530717958648;
w(k) = tmp(k);

for k = 1:36
   % Converge with Halley's iterations, about 5 iterations satisfies
   % the tolerance for most z
   c1 = exp(w);
   c2 = w.*c1 - z;
   w1 = w + (w ~= -1);
   dw = c2./(c1.*w1 - ((w + 2).*c2./(2*w1)));
   w = w - dw;

   if all(abs(dw) < 0.7e-16*(2+abs(w)))
      break;
   end
end

% Specially handle z = 0
w(b ~= 0 & z == 0) = -inf;

% function w = lambertw_aps(branch,x)
% % Lambert_W  Functional inverse of x = w*exp(w).
% % w = Lambert_W(x), same as Lambert_W(x,0)
% % w = Lambert_W(x,0)  Primary or upper branch, W_0(x)
% % w = Lambert_W(x,-1)  Lower branch, W_{-1}(x)
% %
% % See: https://blogs.mathworks.com/cleve/2013/09/02/the-lambert-w-function/
% 
% % Copyright 2013 The MathWorks, Inc.
% 
% % Effective starting guess
% if nargin < 2 || branch ~= -1
%     x = branch;
%     w = ones(size(x));  % Start above -1
% else  
%     w = -2*ones(size(x));  % Start below -1
% end
% v = inf*w;
% 
% % Haley's method
% while any(abs(w - v)./abs(w) > 1.e-8)
%    v = w;
%    e = exp(w);
%    f = w.*e - x;  % Iterate to make this quantity zero
%    w = w - f./((e.*(w+1) - (w+2).*f./(2*w+2)));
% end

% Copyright (C) 1998 Nicol N. Schraudolph <schraudo@inf.ethz.ch>
% Copyright (C) 2016 Colin B. Macdonald
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.

% -*- texinfo -*-
% @documentencoding UTF-8
% @defun  lambertw (@var{z})
% @defunx lambertw (@var{n}, @var{z})
% Compute the Lambert W function of @var{z}.
%
% This function satisfies W(z).*exp(W(z)) = z, and can thus be used to express
% solutions of transcendental equations involving exponentials or logarithms.
%
% @var{n} must be integer, and specifies the branch of W to be computed;
% W(z) is a shorthand for W(0,z), the principal branch.  Branches
% 0 and -1 are the only ones that can take on non-complex values.
%
% For example, the principal branch passes through the point (0, 0):
% @example
% @group
% lambertw (0)
%   @result{} ans = 0
% @end group
% @end example
% And the 0 and -1 branches coincide for the following real value:
% @example
% @group
% x = -1/exp (1);
% lambertw (x)
%   @result{} ans = -1
% lambertw (-1, x)
%   @result{} ans = -1
% @end group
% @end example
%
% If either @var{n} or @var{z} are non-scalar, the function is mapped to each
% element; both may be non-scalar provided their dimensions agree.
% For example, we can repeat the above calculation as:
% @example
% @group
% lambertw ([0 -1], x)
%   @result{} ans =
%       -1  -1
% @end group
% @end example
%
% This implementation should return values within 2.5*eps of its
% counterpart in Maple V, release 3 or later.  Please report any
% discrepancies to the author, Nici Schraudolph <schraudo@@inf.ethz.ch>.
%
% For further algorithmic details, see:
%
% Corless, Gonnet, Hare, Jeffrey, and Knuth (1996), `On the Lambert
% W Function', Advances in Computational Mathematics 5(4):329-359.
%
% @seealso{@@sym/lambertw}
% @end defun

% function w = lambertw(b,z)
%     if (nargin == 1)
%         z = b;
%         b = 0;
%     else
%         % some error checking
%         if (nargin ~= 2)
%             print_usage ();
%         else
%             if (any(round(real(b)) ~= b))
%                 error('branch number for lambertw must be integer')
%             end
%         end
%     end
% 
%     % series expansion about -1/e
%     %
%     % p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
%     % w = (11/72)*p;
%     % w = (w - 1/3).*p;
%     % w = (w + 1).*p - 1
%     %
%     % first-order version suffices:
%     %
%     w = (1 - 2*abs(b)).*sqrt(2*exp(1)*z + 2) - 1;
% 
%     % asymptotic expansion at 0 and Inf
%     %
%     v = log(z + (z == 0 | b == 0)) + 2*pi*1i*b;
%     v = v - log(v + (v == 0));
% 
%     % choose strategy for initial guess
%     %
%     c = abs(z + 1/exp(1));
%     c = (c > 1.45 - 1.1*abs(b));
%     c = c | (b.*imag(z) > 0) | (~imag(z) & (b == 1));
%     w = (1 - c).*w + c.*v;
% 
%     % Halley iteration
%     %
%     for n = 1:10
%         p = exp(w);
%         t = w.*p - z;
%         f = (w ~= -1);
%         t = f.*t./(p.*(w + f) - 0.5*(w + 2.0).*t./(w + f));
%         w = w - t;
%         done = (abs (real (t)) < (2.48*eps)*(1.0 + abs (real (w))) & ...
%                 abs (imag (t)) < (2.48*eps)*(1.0 + abs (imag (w))));
%         if (all (done))
%           break
%         end
%     end
% 
%     % special treatment for infinity and nan
%     isnaninf = (isinf (z) & isreal (z)) | isnan (z);
% 
%     % don't show the warning if we're going to overwrite that entry
%     if (~ all (done | isnaninf))
%       warning ('iteration limit reached, result of lambertw may be inaccurate');
%     end
% 
%     if (any (isnaninf))
%       % broadcast z and b to same size
%       if (isscalar (z) && ~isscalar (b))
%         z = repmat (z, size (b));
%       elseif (~isscalar (z) && isscalar (b))
%         b = repmat (b, size (z));
%       end
% 
%       w(isnan (z)) = nan;
% 
%       I = isinf (z) & isreal (z) & (z > 0);
%       w(I) = inf + 2*b(I)*pi*1i;
% 
%       I = isinf (z) & isreal (z) & (z < 0);
%       w(I) = inf + (2*b(I) + 1)*pi*1i;
% 
%     end
% end


%!assert (isequal (lambertw (0), 0))
%!assert (isequal (lambertw (0, 0), 0))

%!assert (lambertw (-1/exp(1)), -1, 2*eps)
%!assert (lambertw (0, -1/exp(1)), -1, 2*eps)
%!assert (lambertw (-1, -1/exp(1)), -1, 2*eps)

%!test
%! x = [1 2 3 pi 10 100 1000 12345];
%! W = lambertw (x);
%! assert (W.*exp (W), x, -3*eps)

%!test
%! x = [1 2 3 pi 10 100 1000 12345];
%! k = [-3 -2 -1 0 1 2 3 4];
%! W = lambertw (k, x);
%! assert (W.*exp (W), x, -10*eps)

%!test
%! % input shape preserved
%! x = [0 1; 2 3];
%! b = x;
%! W = lambertw (b, x);
%! assert (W.*exp (W), x, -10*eps)

%!test
%! % input shape preserved
%! x = [0 1; 2 3];
%! b = 0;
%! W = lambertw (b, x);
%! assert (W.*exp (W), x, -10*eps)

%!test
%! % input shape preserved
%! x = 10;
%! b = [0 1; 2 3];
%! W = lambertw (b, x);
%! assert (W.*exp (W), x*ones (size (b)), -10*eps)

%!assert (isnan (lambertw (nan)))

%!test
%! % limiting behaviour as z large
%! k = 3;
%! A = lambertw (k, 1e100);
%! assert (abs (imag (A) - 2*pi*k) < 0.1)

%!test
%! % limiting behaviour as z large, up imag axis
%! k = 1;
%! A = lambertw (k, 1e100*1i);
%! assert (abs (imag (A) - (2*k+0.5)*pi) < 0.1)

%!test
%! % limiting behaviour as z large, down imag axis
%! k = -2;
%! A = lambertw (k, -1e100*1i);
%! assert (abs (imag (A) - (2*k-0.5)*pi) < 0.1)

%!test
%! % limiting behaviour as z large, near branch
%! k = 3;
%! A = lambertw (k, -1e100);
%! B = lambertw (k, -1e100 + 1i);
%! C = lambertw (k, -1e100 - 1i);
%! assert (abs (imag (A) - (2*k+1)*pi) < 0.1)
%! assert (abs (imag (B) - (2*k+1)*pi) < 0.1)
%! assert (abs (imag (C) - (2*k-1)*pi) < 0.1)

%!test
%! % infinities and nan
%! A = lambertw ([inf exp(1) -inf nan]);
%! B = [inf  1  inf + pi*1i nan];
%! assert (isequaln (A, B))

%!test
%! % infinities and nan
%! A = lambertw (3, [inf 1 -inf nan]);
%! B = [inf + 2*3*pi*1i  lambertw(3,1)  inf + (2*3+1)*pi*1i  nan];
%! assert (isequaln (A, B))

%!test
%! % infinities and nan
%! A = lambertw ([0 1 2 0], [inf -inf nan exp(1)]);
%! B = [inf  inf+3*pi*1i  nan  1];
%! assert (isequaln (A, B))

%!test
%! % scalar infinity z, vector b
%! A = lambertw ([1 2 -3], inf);
%! B = [lambertw(1, inf)  lambertw(2, inf)  lambertw(-3, inf)];
%! assert (isequal (A, B))

%!test
%! % scalar -infinity z, vector b
%! A = lambertw ([1 2 -3], -inf);
%! B = [lambertw(1, -inf)  lambertw(2, -inf)  lambertw(-3, -inf)];
%! assert (isequal (A, B))

%!test
%! % scalar z nan, vector b
%! A = lambertw ([1 2 -3], nan);
%! B = [nan nan nan];
%! assert (isequaln (A, B))

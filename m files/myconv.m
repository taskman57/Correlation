function c = myconv(a,b)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This function is a low level
  % implementation of conv function
  % to show lower level in FPGA
  % implementation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  a_was_row = isrow(a);
  a = a(:).';                                           % force row
  b = b(:).';                                           % force row
  lena = length(a);
  lenb = length(b);
  conv_len = lenb + lena-1;
  c = zeros(1, conv_len);
  if conv_len -1 > 0                                    % a & b are vectors
      for i=1:conv_len                                  % number of shifts
        min_lim = max(1,(i+1-lenb));
        max_lim = min(i, lena);
        for j=min_lim:max_lim                           % MAC part
          c(i) += a(j)*b(i-j+1);
        endfor
      endfor
  else
    c= b*a;
  endif
  if ~a_was_row
    c = c.';                                            % return column if
  endif                                                 % input was column
endfunction

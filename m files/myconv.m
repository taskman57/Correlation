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
  conv_len = length(b) + length(a)-1;
  zero_pad_a = length(a) - 1;
  zero_pad_b = length(b) - 1;
  printf(' flipped and zero padded.\n\r');
  a = [zeros(1,zero_pad_b) a];
  b = [fliplr(b) zeros(1,zero_pad_a)];
  c = zeros(1, conv_len);
  if conv_len -1 > 0                                    % a & b are vectors
      for i=1:conv_len                                  % number of shifts
        for j=1:length(b)                               % MAC part
          c(i) = c(i) + a(j)*b(j);
        endfor
        b = [0 b(1:end-1)];                             % shift b to right
      endfor
  else
    c= b*a;
  endif
  if ~a_was_row
    c = c.';                                            % return column if
  endif                                                 % input was column
endfunction

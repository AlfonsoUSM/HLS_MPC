%% Write float litarals for C++
% ===============================================================================
% Alfonso Cortes Neira - Universidad Técnica Federico Santa María
% 25-07-2023
% Prints single type values in float format to initialize C++ floats
% ===============================================================================

function str = cpp_float(num)
    if (num == 0)
        str = '     0x0.0p0';
    else
        bin_num = typecast(single(num), 'uint32');
        bin_fra = mod(bin_num, 8388608);
        fraction = dec2hex(2*bin_fra + 16777216);
        bin_exp = idivide(mod(bin_num,2147483648),8388608);
        exponent = int2str(cast(bin_exp, 'int32')-127);
        %if (bin_exp < 127)
        %    exponent = ['-', dec2hex(127 - bin_exp)];
        %else
        %    exponent = [dec2hex(bin_exp - 127), ' '];
        %end
        if (sign(num) == 0)
            str = ['-', '0x1.', fraction(2:7), 'p', exponent];
        else
            str = ['0x1.', fraction(2:7), 'p', exponent];
        end
    end
end
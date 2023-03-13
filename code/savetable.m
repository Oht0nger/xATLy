function savetable(pcdata,outfile)
%SAVE2TXT Save data table as ASCII (.csv or .txt)
%
%   Input
%       pcdata - table object
%       outfile - full output filename

%  Output
%      textfile written to outFile
%
% Lonesome Malambo 08/8/2020, Texas A&M Univeristy

% write table
try
    writetable(pcdata, outfile);
catch ME
    disp(ME.message)
end

end


function installHelpersOnWindows()

% for installing and compiling the helper c-functions (ONLY). For
% the fast multithreading code use LINUX or (try to) compile with opencv
% 
% For this script to work install a compiler, e.g. MINGW see  http://www.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html
% and set the environment variable MW_MINGW64_LOC correctly, see 
% http://www.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html
%
%
cd src

if ~ispc()
    error('This file is only intended for Windows. Use ''make'' for other systems')
end
fnames = {};
targets = {};
HELPER = '..\+fish\+helper\';

fnames{end+1} = 'backtrace_.c';
targets{end+1} = '..\+fish\+core\@DAGraph\private\';

fnames{end+1} = 'getCurrentTracks_.c';
targets{end+1} = '..\+fish\@Tracker\private\';

fnames{end+1} = 'pdist2CenterLine.c';
targets{end+1} = HELPER;

fnames{end+1} = 'pdist2Euclidean.c';
targets{end+1} = HELPER;

fnames{end+1} = 'strucarr2strucmat.c';
targets{end+1} = HELPER;

for i = 1:length(fnames)
   mex('-largeArrayDims',fnames{i},'-outdir',targets{i});
end

% compile CXX
mex -largeArrayDims -cxx CXXFLAGS='-std=c++11 -lstdc++' assignDetectionsToTracks.cpp -outdir ..\+fish\+helper\

cd ..
fprintf('FINISHED!\n Run >> xy.Tracker.runTest() to test the tracking system.\n\n')




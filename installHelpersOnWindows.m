function installHelpersOnWindows(opencvinstallpath)

% for installing and compiling the helper c-functions (ONLY). For
% the fast multithreading code use LINUX or (try to) compile with opencv
% 
% For this script to work install a compiler, e.g. MINGW see  http://www.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html
% and set the environment variable MW_MINGW64_LOC correctly, see 
% http://www.mathworks.com/help/matlab/matlab_external/compiling-c-mex-files-with-mingw.html
%
%
cd src
if nargin<1 || ~exist(opencvinstallpath)
    opencvinstallpath = '';
end

% check for mexopencv
if ~isempty(opencvinstallpath) 
    try 
      evalc('cv.getBuildInformation');
    catch
      error('mexopencv not in path')
    end
    
    
    [a] = fileparts(which('cv.getBuildInformation'));
    mexopencvincludepath = [a '\..\include'];
    mexopencvlibpath = [a '\..\lib']; 
    % find library
    if strfind(computer('arch'),'64')
        opencvlibpath = [opencvinstallpath '\x64'];
    else
        opencvlibpath = [opencvinstallpath '\x86'];
    end
 
    % next compiler name
    d = dir(opencvlibpath);
    idx = [];
    for i = 1:length(d)
        if d(i).isdir && d(i).name(1)~= '.'
           idx = [idx i];
        end
    end
    
    if length(idx)==1
        opencvlibpath = [opencvlibpath '\' d(idx).name '\lib'];
    else 
       error('Cannot establish opencv libpath');
    end
    
    dlibs = dir([opencvlibpath '\opencv*.lib']);
    
   
    opencvincludepath = [opencvinstallpath '\include'];
    % compile into video capture path
    mexargs = {};
    mexargs{end+1} = '-largeArrayDims';
    mexargs{end+1} = '-cxx';
    mexargs{end+1} = 'CXXFLAGS=''-std=c++11 -lstdc++ -fPIC''';
    mexargs{end+1} = ['-I' mexopencvincludepath];
    mexargs{end+1} = ['-I' opencvincludepath];
    mexargs{end+1} = ['-L', opencvlibpath];
    for i = 1:length(dlibs)
        [~,b,~] = fileparts(dlibs(i).name);
        mexargs{end+1} = ['-l' b];
    end
    mexargs{end+1} =  [mexopencvlibpath '\MxArray.obj'];
    
    % videoCapture
    mex(mexargs{:}, 'xyVideoCapture_.cpp',...
        '-outdir', '..\+xy\+core\@VideoCapture\private');
    % VideoHandler
    mex('-c',mexargs{:}, 'VideoHandler.cpp',...
        '-outdir', '..\+xy\+core\@VideoHandlerMex\private');
     
    % videoReader
    mexargs{end+1} = '..\+xy\+core\@VideoHandlerMex\private\VideoHandler.obj';
    mex(mexargs{:}, 'xyVideoHandler_.cpp',...
        '-outdir', '..\+xy\+core\@VideoHandlerMex\private');
   
else
  error('Please make sure that mexopencv toolbox is installed and in the path')
end

if ~ispc()
    error('This file is only intended for Windows. Use ''make'' for other systems')
end
fnames = {};
targets = {};
HELPER = '..\+xy\+helper\';

fnames{end+1} = 'backtrace_.c';
targets{end+1} = '..\+xy\+core\@DAGraph\private\';

fnames{end+1} = 'getCurrentTracks_.c';
targets{end+1} = '..\+xy\@Tracker\private\';

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
mex -largeArrayDims -cxx CXXFLAGS='-std=c++11 -lstdc++' assignDetectionsToTracks.cpp -outdir ..\+xy\+helper\

cd ..
fprintf('FINISHED!\n Run >> xy.Tracker.runTest() to test the tracking system.\n\n')




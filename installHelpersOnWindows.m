function installHelpersOnWindows(opencvflag,opencvinstallpath,flycaptureflag,flycapturelibpath)
% for installing and compiling the c-functions on windows systems. 
% Expects working openCV and mexopencv installation. Please
% consult readme.md for installation instructions. Does not compile the
% flyCapture SDK.

cd src
if nargin<1
    opencvflag = true;
end

if nargin<2 || ~exist(opencvinstallpath)
    opencvinstallpath = '';
end

if nargin<3 
    flycaptureflag = true;
end

if nargin<4 || ~exist(flycapturelibpath)
    flycapturelibpath = '';
end


% check for mexopencv
if opencvflag 
    try 
      evalc('cv.getBuildInformation');
    catch
      error('mexopencv not in path')
    end
    
    if isempty(opencvinstallpath)
        f = regexp(cv.getBuildInformation,'Install path:\s*(?<path>[^\s]+\>)','names');
        if isempty(f)
            error('Cannot determine opencvinstallpath. Please specify');
        else
            opencvinstallpath = f(1).path;
        end
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
    
    fprintf('Using OpenCV Libs installed at "%s"\n',opencvlibpath);
    fprintf('Using MexOpenCV Libs installed at "%s"\n',mexopencvlibpath);
    
    % check for FlyCapture
    if flycaptureflag
        if isempty(flycapturelibpath)
            % assume standard location
            flycapturelibpath = 'C:\Program Files\Point Grey Research\FlyCapture2\lib64';
        end
        if ~exist([flycapturelibpath,filesep,'FlyCapture2.lib'])
            fprintf('**Please provide path to "FlyCapture2.lib". Compiling without PtGray camera support...\n');
            flycaptureflag = false;
        else
            fprintf('Using Flycapture libs in "%s"\n',flycapturelibpath);
        end
        
    end
    
    opencvincludepath = [opencvinstallpath '\include'];
    
    
    % compile into video capture path
    mexargs = {};
    mexargs{end+1} = '-largeArrayDims';
    mexargs{end+1} = '-cxx';
    mexargs{end+1} = 'CXXFLAGS=''-std=c++11 -lstdc++ -fPIC''';
    mexargs{end+1} = ['-I"' mexopencvincludepath,'"',];
    mexargs{end+1} = ['-I"' opencvincludepath,'"'];
    mexargs{end+1} = ['-L"', opencvlibpath,'"'];
    if flycaptureflag
        mexargs{end+1} = ['-L',flycapturelibpath];
        mexargs{end+1} = '-lFlyCapture2';
        mexargs{end+1} = ['-I"',flycapturelibpath, '\..\include"'];
        mexargs{end+1} = '-DFLYCAPTURE';
    end
    
    
    for i = 1:length(dlibs)
        [~,b,~] = fileparts(dlibs(i).name);
        mexargs{end+1} = ['-l' b];
    end
    mexargs{end+1} =  [mexopencvlibpath '\MxArray.obj'];
    
    
    % videoCapture
    mex(mexargs{:}, 'xyVideoCapture_.cpp',...
        '-outdir', '..\+xy\+core\@VideoCapture\private');
    
    if flycaptureflag
        % VideoSaverObj
        mex('-c',mexargs{:},'SaveVideoClass.cpp',...
            '-outdir','..\+xy\+core\@VideoHandlerMex\private');
        mexargs{end+1} = '..\+xy\+core\@VideoHandlerMex\private\SaveVideoClass.obj';
    end
    
   
    % VideoHandler
    mex('-c',mexargs{:}, 'VideoHandler.cpp',...
        '-outdir', '..\+xy\+core\@VideoHandlerMex\private');
     
    % videoReader
    mexargs{end+1} = '..\+xy\+core\@VideoHandlerMex\private\VideoHandler.obj';
    mex(mexargs{:}, 'xyVideoHandler_.cpp',...
        '-outdir', '..\+xy\+core\@VideoHandlerMex\private');
   
    %if flycaptureflag
    %    copyfile([flycapturelibpath '\FlyCapture2.lib'],'..\+xy\+core\@VideoHandlerMex\private\FlyCapture2.lib')
    %end
    
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
fprintf('\nFINISHED!\nRun >> xy.Tracker.runSimpleTest() to test the tracking system.\n\n')




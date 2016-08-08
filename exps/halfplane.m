clear all;

vid = '/data/videos/halfplane/halfplane1.avi'
nbody = 2

xyT = xy.Tracker({0,vid},'useKNN',0,'nbody',nbody, ...
                 'detector.inverted',1,'stmif',1,...
                 'stimulusPresenter',xy.stimulus.PresenterPlane,...
                 'bodylength',130,'bodywidth',30,'displayif',1);


xyT.stimulusPresenter.stmInterval = 30;
xyT.stimulusPresenter.switchInterval = 5;
xyT.stimulusPresenter.col2 = [1,1,1];
xyT.stimulusPresenter.adaptationTime = 0;

% start stim
xyT.track();

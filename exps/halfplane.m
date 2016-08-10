clear all;

vid = '/data/videos/halfplane/halfplane1.avi'
nbody = 2

T = xy.Tracker({0,vid},'useKNN',0,'nbody',nbody, ...
                 'detector.inverted',1,'stmif',1,...
                 'stimulusPresenter',xy.stimulus.PresenterPlane,...
                 'bodylength',130,'bodywidth',30,'displayif',1);


T.stimulusPresenter.stmInterval = 30;
T.stimulusPresenter.switchInterval = 5;
T.stimulusPresenter.col2 = [1,1,1];
T.stimulusPresenter.adaptationTime = 0;

% start stim
T.track();

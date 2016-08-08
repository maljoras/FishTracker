clear all;

vid = '/data/videos/halfplane/halfplane1.avi'
nanimals = 2

ft = xy.Tracker({0,vid},'useKNN',0,'nanimals',nanimals, ...
                 'detector.inverted',1,'stmif',1,...
                 'stimulusPresenter',xy.stimulus.PresenterPlane,...
                 'fishlength',130,'fishwidth',30,'displayif',1);


ft.stimulusPresenter.stmInterval = 30;
ft.stimulusPresenter.switchInterval = 5;
ft.stimulusPresenter.col2 = [1,1,1];
ft.stimulusPresenter.adaptationTime = 0;

% start stim
ft.track();

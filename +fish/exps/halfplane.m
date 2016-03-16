clear all;

vid = '/data/videos/halfplane/halfplane1.avi'
nfish = 2

ft = FishTracker({0,vid},'useKNN',0,'nfish',nfish, ...
                 'detector.inverted',1,'stmif',1,...
                 'stimulusPresenter',FishStimulusPresenterPlane,...
                 'fishlength',130,'fishwidth',30,'displayif',1);


ft.stimulusPresenter.stmInterval = 30;
ft.stimulusPresenter.switchInterval = 5;
ft.stimulusPresenter.col2 = [1,1,1];
ft.stimulusPresenter.adaptationTime = 0;

% start stim
ft.track();

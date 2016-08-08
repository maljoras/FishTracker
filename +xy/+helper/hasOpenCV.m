function bool = hasOpenCV()

  bool = ~~exist('cv.VideoCapture');
  bool = bool && xy.core.VideoCapture.installed();
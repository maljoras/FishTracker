function bool = hasOpenCV()

  bool = ~~exist('cv.VideoCapture');
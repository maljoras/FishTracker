function bool = hasOpenCV()

  bool = ~~exist('cv.VideoCapture');
  bool = bool && fish.core.FishVideoCapture.installed();
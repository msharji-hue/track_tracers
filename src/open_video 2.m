function v = open_video(videoPath)
    v = VideoReader(videoPath);
    fprintf('Loaded: %d frames at %.2f FPS\n', floor(v.Duration * v.FrameRate), v.FrameRate);
end

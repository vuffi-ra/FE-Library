function res = run_all()
  [~, folder]    = fileparts(fileparts(mfilename('fullpath')));
  test_files     = dir(fullfile(folder, 'test*'));
  test_files_str = fullfile(folder, {test_files.name});
  res            = runtests(test_files_str);
end

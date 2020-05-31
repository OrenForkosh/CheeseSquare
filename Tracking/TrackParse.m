function data = TrackParse(obj)

data.Type = regexprep(obj.FilePrefix, '^([^\.]*)\..*', '$1');
data.GroupID = str2double(regexprep(obj.FilePrefix, '.*exp([0-9]*).*', '$1'));
data.Day = str2double(regexprep(obj.FilePrefix, '.*day([0-9]*).*', '$1'));
data.CameraID = str2double(regexprep(obj.FilePrefix, '.*cam([0-9]*).*', '$1'));

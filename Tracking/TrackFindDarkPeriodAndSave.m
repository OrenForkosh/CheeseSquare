function obj=TrackFindDarkPeriodAndSave(obj)
% Auxilary function. Finds dark\light transitions (using
% TrackFindDarkPeriod) and saves the data back to 'obj'
obj = TrackLoad(obj);
obj = TrackFindDarkPeriod(obj, .2);
TrackSave(obj);
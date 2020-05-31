function obj = TrackExport(obj)
% Export only the tracking data to file. This data is primarily used for 
% displaying the videos with overlaid trajectories (from BigCheese)
obj = TrackLoad(obj);
%%
export.x = double(obj.x);
export.y = double(obj.y);
export.time = obj.time;
export.colors = double(obj.Colors.Centers);
export.meta.subject.centerColors = double(obj.Colors.Centers);
export.zones = double(obj.zones);
export.labels = obj.ROI.ZoneNames;
export.nSubjects = obj.nSubjects;
%%
% export.Events(1).Title = 'contacts';
% idx = 1;
% events.Begin = [];
% events.End = [];
% events.Desc = [];
% events.Members = [];
% for s1=1:obj.nSubjects
%     for s2=s1+1:obj.nSubjects
%         for i=1:length(obj.Hierarchy.ChaseEscape.interactions.events{s1, s2})
%             export.Events(1).instances(idx).Begin = obj.Hierarchy.ChaseEscape.interactions.events{s1, s2}{i}.BegFrame / obj.FrameRate;
%             export.Events(1).instances(idx).End = obj.Hierarchy.ChaseEscape.interactions.events{s1, s2}{i}.EndFrame / obj.FrameRate;
%             export.Events(1).instances(idx).Desc = obj.Hierarchy.ChaseEscape.interactions.events{s1, s2}{i}.desc;
%             export.Events(1).instances(idx).Members = uint8([s1, s2]);
%         end
%     end
% end

%%
old = true;
if false;
try
    if old
        export.Events(1).Title = 'contacts';
        idx = 1;
        events.Begin = [];
        events.End = [];
        %events.Desc = {};
        events.Members = [];
        for s1=1:obj.nSubjects
            for s2=s1+1:obj.nSubjects
                for i=1:length(obj.Hierarchy.ChaseEscape.interactions.events{s1, s2})
                    events.Begin(idx) = obj.Hierarchy.ChaseEscape.interactions.events{s1, s2}{i}.BegFrame / obj.FrameRate;
                    events.End(idx) = obj.Hierarchy.ChaseEscape.interactions.events{s1, s2}{i}.EndFrame / obj.FrameRate;
                    %            events.Desc{idx} = obj.Hierarchy.ChaseEscape.interactions.events{s1, s2}{i}.desc;
                    %events.Desc{idx} = '';
                    events.Members(idx) = s1 * 10 + s2;
                    idx = idx + 1;
                end
            end
        end
        [qqq, o] = sort(events.End);
        events.Begin = events.Begin(o);
        events.End = events.End(o);
        %events.Desc = {events.Desc{o}};
        events.Members = events.Members(o);
        
        % %export.Events(1).instances = struct('Begin', [], 'End', [], 'Desc', [], 'Members', []);
        % export.Events(1).instances = struct('Begin', [], 'End', [], 'Members', []);
        % for i=1:min(length(events.Begin), 500)
        %     export.Events(1).instances(i).Begin = events.Begin(i);
        %     export.Events(1).instances(i).End = events.End(i);
        %     %export.Events(1).instances(i).Desc = events.Desc{i};
        %     export.Events(1).instances(i).Members = events.Members(i);
        % end
        export.Events(1).data = [events.Begin; events.End; events.Members];
    else
        export.Events(1).Title = 'contacts';
        events.Begin = min([obj.Contacts.List.beg]);
        events.End = max([obj.Contacts.List.end]);
        events.Members = zeros(1, length(obj.Contacts.List));
        for i=1:length(obj.Contacts.List)
            events.Members(i) = 10.^[0:length(obj.Contacts.List(i).subjects)-1] * obj.Contacts.List(i).subjects';
        end
        export.Events(1).data = [events.Begin; events.End; events.Members];
    end
end
end
%%
filename = [regexprep(obj.VideoFile, '\.[^\.]*$', '') '.mat'];
social = export;
save(filename, 'social', '-v6');


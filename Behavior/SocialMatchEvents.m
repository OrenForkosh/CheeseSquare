function obj = SocialMatchEvents(obj, events)
% Match interaction events between subjects (as interactions are computed
% separetly for each mouse)
%
%       Created by OREN FORKOSH
%

newe = [];
curr = struct();
idx = 1;
for s1=1:obj.nSubjects    
    for s2=s1+1:obj.nSubjects
        if s1==s2; continue; end
        e1 = events{s1, s2};
        e2 = events{s2, s1};
        e = [e2{:}];
        if ~isfield(e, 'EndFrame')
            continue;
        end
        endf = [e.EndFrame];
        begf = [e.BegFrame];
        %%
%         
%         bf1 = [E1.BegFrame];
%         ef1 = [E1.EndFrame];
% 
%         E2 = [e2{:}];
%         bf2 = [E2.BegFrame];
%         ef2 = [E2.EndFrame];
        
        %%
        E1 = [e1{:}];
        E2 = [e2{:}];
        overlap = Q.overlap([E1.BegFrame], [E1.EndFrame], [E2.BegFrame], [E2.EndFrame]);
        for i=1:size(overlap, 1)
            for j=find(overlap(i, :))
                curr.beg = [e1{i}.BegFrame; begf(j)];
                curr.end = [e1{i}.EndFrame; endf(j)];
                curr.valid = true; %e1{i}.valid || e(j).valid;
                curr.data{1} = e1{i}.data;
                curr.data{2} = e(j).data;
                curr.subjects = [s1 s2]';
                curr.states = false(2, length(obj.Interactions.PredPrey.model.names));
                curr.states(1, unique(curr.data{1})) = true;
                curr.states(2, unique(curr.data{2})) = true;
                
                if isempty(newe)
                    newe = curr;
                else
                    newe(idx) = curr;
                end
                idx = idx + 1;
            end
        end
        %%
%         %subject = [ones(1, length(e1)) * s1 ones(1, length(e2)) * s2];
%         for i=1:length(e1)
%             m = ~(e1{i}.BegFrame > endf | e1{i}.EndFrame < begf);
%             for j=find(m)
%                 curr.beg = [e1{i}.BegFrame; begf(j)];
%                 curr.end = [e1{i}.EndFrame; endf(j)];
%                 curr.valid = true; %e1{i}.valid || e(j).valid;
%                 curr.data{1} = e1{i}.data;
%                 curr.data{2} = e(j).data;
%                 curr.subjects = [s1 s2]';
%                 curr.states = false(2, length(obj.Interactions.PredPrey.model.names));
%                 curr.states(1, unique(curr.data{1})) = true;
%                 curr.states(2, unique(curr.data{2})) = true;
%                 
%                 if isempty(newe)
%                     newe = curr;
%                 else
%                     newe(idx) = curr;
%                 end
%                 idx = idx + 1;
%             end
%         end
    end
end
[v, o] = sort(min([newe.beg]));
obj.Contacts.List = newe(o);
valid = [obj.Contacts.List.valid];
obj.Contacts.List = obj.Contacts.List(valid);
obj.Contacts.StateNames = obj.Interactions.PredPrey.model.names;
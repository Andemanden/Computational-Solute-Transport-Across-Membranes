classdef Kontrolvol %Controlvolume class for intuity and e
   properties
      Prop=0;
      timeinterval=0;
   end
   methods
       function obj = KontrolClass(val)
         if nargin > 0
            obj.Prop = val;
            %disp(obj.Prop);
            obj.timeinterval=val;
         end
      end
   end
end




classdef Kontrolvol
   properties
      Prop=0;
   end
   methods
      function obj = MyClass(val)
         if nargin > 0
            obj.Prop = val;
            disp(obj.Prop);
         end
      end
   end
end




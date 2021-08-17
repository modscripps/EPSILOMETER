function [trueOrFalse] = isclassfield(classname,fieldname)

trueOrFalse = any(cell2mat(strfind(fields(classname),fieldname)));
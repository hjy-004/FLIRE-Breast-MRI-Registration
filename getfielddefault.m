function val = getfielddefault(s,tag,defval)

if isfield(s,tag)
  val = getfield(s,tag);
else
  val = defval;
end


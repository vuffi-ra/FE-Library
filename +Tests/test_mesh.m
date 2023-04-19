function tests = test_mesh
  tests = functiontests(localfunctions);
end


function test_unitsquare(~)
  mesh = Mesh.UnitSquare();
  mesh.plot();
  mesh.plot(true);
end

function test_refine(~)
  mesh = Mesh.UnitSquare();
  mesh = Mesh.Refine(mesh);
  mesh.plot(true);
end

function test_golfofmexico(~)
  Mesh.GolfOfMexico();
end

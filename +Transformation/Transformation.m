classdef Transformation < handle

  properties (SetAccess = protected, GetAccess = public)
    mesh (1, 1) Mesh.General = Mesh.Empty()
  end

  methods (Access = public)
      function obj = Transformation(mesh_in)
      arguments
        mesh_in (1, 1) Mesh.General
      end

      obj.mesh = mesh_in;
    end

    % Coordinates (xi, eta) in reference triangle
    % Number of element e

    function [x, y] = referenceToMesh(xi, eta, e)
      a = obj.mesh.XY(obj.mesh.V0T(e, 1), :);
      b = obj.mesh.XY(obj.mesh.V0T(e, 2), :);
      c = obj.mesh.XY(obj.mesh.V0T(e, 3), :);

      [x, y] = xi * (a - c) + eta * (b - c) + c;
    end

    function J = referenceToMeshJacobian(~, ~, e)
      a = obj.mesh.XY(obj.mesh.V0T(e, 1), :);
      b = obj.mesh.XY(obj.mesh.V0T(e, 2), :);
      c = obj.mesh.XY(obj.mesh.V0T(e, 3), :);
      J = [a(1) - c(1), b(1) - c(1);
           a(2) - c(2), b(2) - c(2)];
    end
  end
end

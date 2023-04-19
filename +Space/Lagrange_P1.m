classdef Lagrange_P1

  properties (SetAccess = protected, GetAccess = public)
    trafo (1, 1) Trafo.Trafo
  end

  methods (Access = public)
    function obj = Lagrange_P1(trafo_in)
      obj.trafo = trafo_in;
    end

    function v = val(xi, eta, i)
      %% ...
    end

    function v = grad(xi, eta, i)
      %% ...
    end

    function i = dofmapper(hat_i, e)
      %% ...
    end
  end
end

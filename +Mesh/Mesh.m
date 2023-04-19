classdef Mesh < handle

  properties (SetAccess = public, GetAccess = public)
    XY    (:, 2) double %% coordinats
    V0T   (:, 3) double %% table of connectivity
    idE0T (:, 3) double %% type of edge
    E0T   (:, 3) double %% global number of edge
  end

  methods (Access = public)
    
    function obj = Mesh(vertices, faces, edgeTypes, edgeNumbers)
      obj.XY = vertices;
      obj.V0T = faces;
      obj.idE0T = edgeTypes;
      obj.E0T = edgeNumbers;
    end

    function plot(mesh, showDetails)
      arguments
        mesh        (1, 1) Mesh.Mesh
        showDetails (1, 1) logical = false
      end

      hold('on')
      daspect([1, 1, 1]) % adjust aspect ration, requires Octave >= 3.8
      textarray = @(x1,x2,s,color) arrayfun(@(a,b,c) ...
                                             text(a,b,int2str(c),'Color',color,'FontSize',16), x1, x2, s);
      numV = size(mesh.XY, 1);
      numT = size(mesh.V0T,1);

      %% Triangle boundaries.
      trisurf(mesh.V0T,mesh.XY(:,1),mesh.XY(:,2),zeros(numV,1), 'facecolor', 'none', 'LineWidth', 2);
      numE = max(max(mesh.E0T));
      coordVx = mesh.XY(:,1);
      coordVy = mesh.XY(:,2);
      baryT = [ sum(reshape(coordVx(mesh.V0T), size(mesh.V0T, 1), []), 2), ...
                sum(reshape(coordVy(mesh.V0T), size(mesh.V0T, 1), []), 2) ] / 3;
      baryE = zeros(numE,2);
      areaE = zeros(numE,1);
      nuE = zeros(numE,2);
      idE = zeros(numE,1);
      for k = 1:numT
        for i = 1:3
          v1 = mod(i,3)+1;
          v2 = mod(i+1,3)+1;
          edge = mesh.E0T(k,i);
          baryE(edge,:) = 0.5 * (mesh.XY(mesh.V0T(k,v1),:) + mesh.XY(mesh.V0T(k,v2),:) );
          idE(edge) = mesh.idE0T(k,i);
          E = mesh.XY(mesh.V0T(k,v1),:) - mesh.XY(mesh.V0T(k,v2),:);
          areaE(edge) = norm(E);
          nuE(edge,:) = [-E(2) E(1)] / areaE(edge);
        end
      end
      if showDetails
        %% Local edge numbers.
        w = [1/12, 11/24, 11/24; 11/24, 1/12, 11/24; 11/24, 11/24, 1/12];
        for kE = 1 : 3
          textarray(reshape(mesh.XY(mesh.V0T,1),numT,3)*w(:,kE), ...
                    reshape(mesh.XY(mesh.V0T,2),numT,3)*w(:,kE), kE*ones(numT, 1),'blue');
        end % for
        %% Global vertex numbers.
        textarray(mesh.XY(:,1), mesh.XY(:,2), (1:numV)','magenta');
        %% Local vertex numbers.
        w = [5/6, 1/12, 1/12; 1/12, 5/6, 1/12; 1/12, 1/12, 5/6];
        for kV = 1 : 3
          textarray(reshape(mesh.XY(mesh.V0T,1),numT,3)*w(:,kV), ...
                    reshape(mesh.XY(mesh.V0T,2),numT,3)*w(:,kV), kV*ones(numT, 1),'red');
        end % for
        %% Global edge numbers.
        textarray(baryE(:,1), baryE(:,2), (1:numE)','cyan');
        %% Triangle numbers.
        textarray(baryT(:,1), baryT(:,2), (1:numT)','black');
      end % if
      %% Edge IDs.
      markEext = idE ~= 0; % mark boundary edges
      textarray(baryE(markEext,1) + nuE(markEext,1).*areaE(markEext)/8, ...
                baryE(markEext,2) + nuE(markEext,2).*areaE(markEext)/8, idE(markEext),'green');
      minCoordV = min(mesh.XY, [], 1);
      maxCoordV = max(mesh.XY, [], 1);
      dist = maxCoordV - minCoordV;
      space = 0.05*max(dist);
      axis([min(mesh.XY(:,1)) max(mesh.XY(:,1)) min(mesh.XY(:,2)) max(mesh.XY(:,2))]+space*[-1 1 -1 1]);
    end

    function refined = refine(mesh)
      arguments
        mesh (1, 1) Mesh.General
      end

      numV = size(mesh.XY,1);
      numT = size(mesh.V0T,1);
      numE = max(max(mesh.E0T));

      %% ctrV = numV; % vertex counter, set to number of vertices in the coarse mesh
      numT2 = 4*numT;
      refinedXY = [mesh.XY; zeros(numE,2)];

      newVertices   = zeros(numE, 1);
      refinedV0T       = zeros(numT2,3);
      refinedIdE0T     = zeros(numT2,3);
      edgeMidpoints = zeros(1,3);

      for k = 1:numT
        element = 4*(k-1);
        for i = 1:3
          i1       = mod(i  ,3) + 1;
          i2       = mod(i+1,3) + 1;
          midpoint = 0.5 * (mesh.XY(mesh.V0T(k,i1),:) + mesh.XY(mesh.V0T(k,i2),:));
          edge     = mesh.E0T(k,i);
          if newVertices(edge) == 0
            %% ctrV                        = ctrV + 1;
            %% newVertices(edge)           = ctrV;
            newVertices(edge)           = edge + numV;
            refinedXY(newVertices(edge),:) = midpoint;
            if mesh.idE0T(k,i) ~= 0
              refinedIdE0T(element+i1,i) = mesh.idE0T(k,i);
              refinedIdE0T(element+i2,i) = mesh.idE0T(k,i);
            end
          end
          edgeMidpoints(i) = newVertices(edge);
        end
        obj.V0T(element+1,:) = [mesh.V0T(k,1),    edgeMidpoints(3), edgeMidpoints(2)];
        obj.V0T(element+2,:) = [edgeMidpoints(3), mesh.V0T(k,2),    edgeMidpoints(1)];
        obj.V0T(element+3,:) = [edgeMidpoints(2), edgeMidpoints(1), mesh.V0T(k,3)];
        obj.V0T(element+4,:) = [edgeMidpoints(1), edgeMidpoints(2), edgeMidpoints(3)];
      end

      numT = size(refinedV0T,1);
      numV = size(refinedXY,1);
      V2T  = sparse(refinedV0T(:, [1 2 3 1 2 3 1 2 3]), refinedV0T(:, [2 3 1 2 3 1 2 3 1]), ...
                    [(1:numT)', zeros(numT,3),(1:numT)', zeros(numT,3),(1:numT)'], numV, numV);
      [r, c] = find(triu(V2T + V2T'));
      V2E = sparse(r, c, 1:size(r, 1), numV, numV);
      V2E = V2E + V2E';
      refinedE0T = full(V2E(sub2ind([numV,numV], refinedV0T(:,[2,3,1]), refinedV0T(:,[3,1,2]))));

      refined = Mesh.Mesh(refinedXY, refinedV0T, refinedIdE0T, refinedE0T);
    end
  end
end

function buildvoxels
    filename = '04.obj';
    [vertices, faces] = readObj(filename);
    
    voxelSize = 0.1;
    [voxels, voxelGrid] = fitRoomByVoxels(vertices, faces, voxelSize);
    
    figure;
    hold on;
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Voxelized Room');
    
    for i = 1:size(voxels, 1)
        voxel = voxels(i, :);
        drawVoxel(voxel, voxelSize);
    end
end

function [vertices, faces] = readObj(filename)
    fid = fopen(filename, 'r');
    
    vertices = [];
    faces = [];
    
    while ~feof(fid)
        line = fgetl(fid);
        
        if startsWith(line, 'v ')
            vertex = sscanf(line, 'v %f %f %f');
            vertices = [vertices; vertex'];
        elseif startsWith(line, 'f ')
            face = sscanf(line, 'f %d %d %d');
            faces = [faces; face'];
        end
    end
    
    fclose(fid);
end

function [voxels, voxelGrid] = fitRoomByVoxels(vertices, faces, voxelSize)
    minX = min(vertices(:, 1));
    minY = min(vertices(:, 2));
    minZ = min(vertices(:, 3));
    maxX = max(vertices(:, 1));
    maxY = max(vertices(:, 2));
    maxZ = max(vertices(:, 3));
    
    gridSizeX = ceil((maxX - minX) / voxelSize);
    gridSizeY = ceil((maxY - minY) / voxelSize);
    gridSizeZ = ceil((maxZ - minZ) / voxelSize);
    
    voxelGrid = false(gridSizeX, gridSizeY, gridSizeZ);
    
    for i = 1:gridSizeX
        for j = 1:gridSizeY
            for k = 1:gridSizeZ
                voxelX = minX + (i - 0.5) * voxelSize;
                voxelY = minY + (j - 0.5) * voxelSize;
                voxelZ = minZ + (k - 0.5) * voxelSize;
                voxelPosition = [voxelX, voxelY, voxelZ];
                
                inside = isPointInsideRoom(voxelPosition, vertices, faces);
                
                voxelGrid(i, j, k) = inside;
            end
        end
    end
    
    [voxelsX, voxelsY, voxelsZ] = ind2sub(size(voxelGrid), find(voxelGrid));
    voxels = [voxelsX, voxelsY, voxelsZ];
end

function inside = isPointInsideRoom(point, vertices, faces)
    intersectionCount = 0;
    
    for i = 1:size(faces, 1)
        face = faces(i, :);
        
        vertex1 = vertices(face(1), :);
        vertex2 = vertices(face(2), :);
        vertex3 = vertices(face(3), :);
        
        if rayIntersectsTriangle(point, vertex1, vertex2, vertex3)
            intersectionCount = intersectionCount + 1;
        end
    end
    
    inside = mod(intersectionCount, 2) == 1;
end

function intersects = rayIntersectsTriangle(origin, vertex1, vertex2, vertex3)
    edge1 = vertex2 - vertex1;
    edge2 = vertex3 - vertex1;
    normal = cross(edge1, edge2);
    
    if abs(dot(normal, origin - vertex1)) < eps
        intersects = false;
        return;
    end
    
    t = dot(normal, vertex1 - origin) / dot(normal, origin);
    intersectionPoint = origin + t * (origin - vertex1);
    
    edge1Test = cross(edge1, intersectionPoint - vertex1);
    edge2Test = cross(edge2, intersectionPoint - vertex2);
    intersects = dot(edge1Test, normal) >= 0 && dot(edge2Test, normal) >= 0;
end

function drawVoxel(voxel, voxelSize)
    corners = [
        voxel(1) - voxelSize/2, voxel(2) - voxelSize/2, voxel(3) - voxelSize/2;
        voxel(1) + voxelSize/2, voxel(2) - voxelSize/2, voxel(3) - voxelSize/2;
        voxel(1) + voxelSize/2, voxel(2) + voxelSize/2, voxel(3) - voxelSize/2;
        voxel(1) - voxelSize/2, voxel(2) + voxelSize/2, voxel(3) - voxelSize/2;
        voxel(1) - voxelSize/2, voxel(2) - voxelSize/2, voxel(3) + voxelSize/2;
        voxel(1) + voxelSize/2, voxel(2) - voxelSize/2, voxel(3) + voxelSize/2;
        voxel(1) + voxelSize/2, voxel(2) + voxelSize/2, voxel(3) + voxelSize/2;
        voxel(1) - voxelSize/2, voxel(2) + voxelSize/2, voxel(3) + voxelSize/2;
    ];
    
    patch('Vertices', corners, 'Faces', [
        1 2 3 4;
        2 6 7 3;
        6 5 8 7;
        5 1 4 8;
        1 2 6 5;
        4 3 7 8;
    ], 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

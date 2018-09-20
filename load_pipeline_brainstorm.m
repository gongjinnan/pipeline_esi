function [ sMri, Cortex, sHead ] = load_pipeline_brainstorm(MriFile,FileFormatMri,TransFile,TessFileL,TessFileR,FileFormatTess,nVertices,erodeFactor,fillFactor )

sMri = in_mri(MriFile, FileFormatMri, 0);

fid = fopen(TransFile, 'rt');
strFid = fread(fid, [1 Inf], '*char');
fclose(fid);
% Evaluate the file (.m file syntax)
eval(strFid);

Tbst2ft = [diag([-1, 1, 1] ./ sMri.Voxsize), [size(sMri.Cube,1); 0; 0]; 0 0 0 1];
% Set the MNI=>SCS transformation in the MRI
Tmni = transform.vox07mm2spm * Tbst2ft;
sMri.NCS.R = Tmni(1:3,1:3);
sMri.NCS.T = Tmni(1:3,4);
% MNI coordinates for the AC/PC/IH fiducials
AC = [0,   3,  -4] ./ 1000;
PC = [0, -25,  -2] ./ 1000;
IH = [0, -10,  60] ./ 1000;
Origin = [0, 0, 0];
% Convert: MNI (meters) => MRI (millimeters)
sMri.NCS.AC     = cs_convert(sMri, 'mni', 'mri', AC) .* 1000;
sMri.NCS.PC     = cs_convert(sMri, 'mni', 'mri', PC) .* 1000;
sMri.NCS.IH     = cs_convert(sMri, 'mni', 'mri', IH) .* 1000;
sMri.NCS.Origin = cs_convert(sMri, 'mni', 'mri', Origin) .* 1000;

Tscs = transform.vox07mm2bti * Tbst2ft;
sMri.SCS.R = Tscs(1:3,1:3);
sMri.SCS.T = Tscs(1:3,4);
% Standard positions for the SCS fiducials
NAS = [90,   0, 0] ./ 1000;
LPA = [ 0,  75, 0] ./ 1000;
RPA = [ 0, -75, 0] ./ 1000;
Origin = [0, 0, 0];
% Convert: SCS (meters) => MRI (millimeters)
sMri.SCS.NAS    = cs_convert(sMri, 'scs', 'mri', NAS) .* 1000;
sMri.SCS.LPA    = cs_convert(sMri, 'scs', 'mri', LPA) .* 1000;
sMri.SCS.RPA    = cs_convert(sMri, 'scs', 'mri', RPA) .* 1000;
sMri.SCS.Origin = cs_convert(sMri, 'scs', 'mri', Origin) .* 1000;


TessL = in_tess(TessFileL, FileFormatTess, sMri, []);

NewTessL = db_template('surfacemat');
NewTessL.Comment  = TessL(1).Comment;
NewTessL.Vertices = TessL(1).Vertices;
NewTessL.Faces    = TessL(1).Faces;

TessR = in_tess(TessFileR, FileFormatTess, sMri, []);

 NewTessR = db_template('surfacemat');
NewTessR.Comment  = TessR(1).Comment;
NewTessR.Vertices = TessR(1).Vertices;
NewTessR.Faces    = TessR(1).Faces;

% Merge surfaces
Cortex = my_tess_concatenate([NewTessL, NewTessR], sprintf('cortex_%dV', size(NewTessL.Vertices,1) + size(NewTessR.Vertices,1)), 'Cortex');

 [Cortex] = my_in_tess_bst( Cortex, 1 );
sMri.Histogram = mri_histogram(sMri.Cube);


% Threshold mri to the level estimated in the histogram
headmask = (sMri.Cube > sMri.Histogram.bgLevel);
% Closing all the faces of the cube
headmask(1,:,:)   = 0*headmask(1,:,:);
headmask(end,:,:) = 0*headmask(1,:,:);
headmask(:,1,:)   = 0*headmask(:,1,:);
headmask(:,end,:) = 0*headmask(:,1,:);
headmask(:,:,1)   = 0*headmask(:,:,1);
headmask(:,:,end) = 0*headmask(:,:,1);




% Erode + dilate, to remove small components
if (erodeFactor > 0)
    headmask = headmask & ~mri_dilate(~headmask, erodeFactor);
    headmask = mri_dilate(headmask, erodeFactor);
end

headmask = (mri_fillholes(headmask, 1) & mri_fillholes(headmask, 2) & mri_fillholes(headmask, 3));
[sHead.Faces, sHead.Vertices] = isosurface(headmask);
% Downsample to a maximum number of vertices
maxIsoVert = 60000;
if (length(sHead.Vertices) > maxIsoVert)   
    [sHead.Faces, sHead.Vertices] = reducepatch(sHead.Faces, sHead.Vertices, maxIsoVert./length(sHead.Vertices));   
end

% Remove small objects
[sHead.Vertices, sHead.Faces] = tess_remove_small(sHead.Vertices, sHead.Faces);

% Downsampling isosurface
[sHead.Faces, sHead.Vertices] = reducepatch(sHead.Faces, sHead.Vertices, nVertices./length(sHead.Vertices));

% Convert to millimeters
sHead.Vertices = sHead.Vertices(:,[2,1,3]);
sHead.Faces    = sHead.Faces(:,[2,1,3]);
sHead.Vertices = bst_bsxfun(@times, sHead.Vertices, sMri.Voxsize);
% Convert to SCS
sHead.Vertices = cs_convert(sMri, 'mri', 'scs', sHead.Vertices ./ 1000);

% Reduce the final size of the meshed volume
erodeFinal = 3;
% Fill holes in surface
%if (fillFactor > 0)
    
[sHead.Vertices, sHead.Faces] = tess_fillholes(sMri, sHead.Vertices, sHead.Faces, fillFactor, erodeFinal);
 sHead.Atlas = db_template('Atlas');
 sHead.VertConn = tess_vertconn(sHead.Vertices, sHead.Faces);
 sHead.VertNormals = tess_normals(sHead.Vertices, sHead.Faces, sHead.VertConn);
 sHead.Curvature = single(tess_curvature(sHead.Vertices, sHead.VertConn, sHead.VertNormals, .1));
sHead.SulciMap = tess_sulcimap(sHead);


end


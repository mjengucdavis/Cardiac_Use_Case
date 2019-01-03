
figure;

fprintf ( 1, '\n' );
fprintf ( 1, 'MAKE_AVI_MOVIE_EXMPLE2\n' );
fprintf ( 1, '  Create an AVI animation of 2D graphics,\n' );
fprintf ( 1, '  generating one frame at a time.\n' )
%
%  Set the number of frames.
%
numframes = 1000;
%
vidobj = VideoWriter ( '2D_OUTPUT_movie.avi');
open(vidobj);
colormap(jet);
%  Now generate each frame.
%
for jj=1:numframes
    
    str=sprintf('modelOUTPUT_ap%d.dat', jj);
    v=load(str);
    contourf(v,[-100:1:50],'LineStyle','none');
    caxis([-100 50]);
    colorbar;
    
    frame = getframe ( gca );
    writeVideo(vidobj,frame);
    
end

%
%  Once all the frames have been generated and added, we must close the file.
%
close ( vidobj );

fprintf ( 1, '\n' );
fprintf ( 1, 'MAKE_AVI_MOVIE_EXAMPLE2\n' );
fprintf ( 1, '  Normal end of execution.\n' );
fprintf ( 1, '  The movie file "mov_intensity.avi" has been created.\n' );

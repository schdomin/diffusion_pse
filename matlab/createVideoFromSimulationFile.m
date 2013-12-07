function createVideoFromSimulationFile( )
clear;

%ds filepath
strFilepath = '../bin/simulation.txt';

%ds open the file
fileID = fopen( strFilepath );

%ds get the first line
cCell = textscan( fileID, '%u %u', 1 );

%ds get number of particles N and timesteps T
uNumberOfGridPoints1D = cCell{1};
uNumberOfTimesteps    = cCell{2};

%ds informative
disp( ['Number of Grid Points 1D: ', num2str( uNumberOfGridPoints1D ) ] );
disp( ['    Number of Time Steps: ', num2str( uNumberOfTimesteps ) ] );

%ds data structure
matHeatGrid = zeros( uNumberOfTimesteps, uNumberOfGridPoints1D, uNumberOfGridPoints1D );

disp( [ 'starting data import from: ', strFilepath ] ); 
tic;

%ds for each timestep
for uCurrentTimestep = 1:1:uNumberOfTimesteps
   
    %ds scan all lines
    for uCurrentLine = 1:1:uNumberOfGridPoints1D
        
        %ds set the current line
        cCellLine = textscan( fileID, '%f', uNumberOfGridPoints1D );  
        
        matHeatGrid( uCurrentTimestep, uCurrentLine, : ) = cCellLine{:};
        
    end
end

disp( [ 'finished data import - time: ', num2str( toc ) ] );

%ds prepare video writer
writerObj = VideoWriter('simulation.avi');
writerObj.FrameRate = 25;
open( writerObj );
disp( [ 'timestep: ', num2str( 1 ) ] );

%ds create initial data with the first timestep
matInitial = squeeze( matHeatGrid( 1, :, : ) );

%ds get max/min value for constant color mapping
dMaxValue = max( max( matInitial ) );
dMinValue = abs( min( min( matInitial ) ) );

%ds shift the whole matrix so the minimum value is zero
matHeatGrid = matHeatGrid + dMinValue;

%ds also shift the max value
dMaxValue = dMaxValue + dMinValue;

%ds plot first frame
contourf( matInitial, 100, 'LineColor', 'none' );
caxis([ 0 dMaxValue ]);
set( gca, 'nextplot' ,'replacechildren' );
set( gca, 'color', [0.5, 1, 0.5] );
set( gcf, 'Renderer' ,'zbuffer' );
%set( gcf, 'LineStyle', 'none' );

% %ds for each remaining timestep
for uCurrentTimestep = 2:1:uNumberOfTimesteps
    
    %ds only save every 20th frame instead 500 -> 25 frames per sec
    if mod( uCurrentTimestep, 10 ) == 0
        
        %ds create a figure
        contourf( squeeze( matHeatGrid( uCurrentTimestep, :, : ) ), 100, 'LineColor', 'none'  );
        caxis([ 0 dMaxValue ]);
        frame = getframe( gcf );
        writeVideo( writerObj, frame );

        disp( [ 'timestep: ', num2str( uCurrentTimestep ) ] );

    end
end

%ds write video file
close( writerObj );

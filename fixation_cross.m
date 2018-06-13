function rect = fixation_cross(width,height,centrex,centrey)
% $Id: fixation_cross.m 4 2007-06-23 11:13:19Z damienm $
%Make a small fixation cross to display on screen.
%Width and height are in pixels.
%centrex & centrey give the x & y screen coordinates where you want the cross centred.

	rect = zeros(4,2);
	
	width = width/2;
	height = height/2;
	
	rect(1,1) = -height;
	rect(2,1) = width;
	rect(3,1) = height;
	rect(4,1) = -width;
				
	rect(1,2) = -width;
	rect(2,2) = height;
	rect(3,2) = width;
	rect(4,2) = -height;
	
		
	rect(1:2:4,:) = rect(1:2:4,:) + centrex;
	rect(2:2:4,:) = rect(2:2:4,:) + centrey;
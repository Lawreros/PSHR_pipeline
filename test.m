close ALL;
changeSlider;
disp('hey');

function changeSlider
%   fig2 = uifigure('Position',[125 100 300 300]);
%   someth = uibutton(fig2, 'ButtonPushedFcn', @(src,event)randomplot(3));

  fig = uifigure('Position',[0 100 1000 600]);
  fig2 = figure('Position',[0 100 1000 600]);
  ax2 = axes(fig2);
  %ax=uiaxes(fig, 'Position', [50 400 900 150]);
  s = uislider(fig,'Position',[75 150 150 3]);
  incrementSlider;
  b = uibutton(fig,'Position',[100 50 100 22], ...
    'Text','Increment', ...
    'ButtonPushedFcn',@(src,event)incrementSlider(ax2));
  p =uibutton(fig, 'ButtonPushedFcn', @(src,event)randomplot(1,ax2,fig,3)); 



  function incrementSlider(ax2)
    if s.Value < s.Limits(2)
      s.Value = s.Value + 1;
    end    
  end

    function randomplot(entries,axes,fi,timer)
       while timer < 10
        plot(axes, rand(20,entries));
        %q = uibutton(fi, 'Position', [20,20,20,20], 'Text','hey my dude','ButtonPushedFcn','disp("gotit")');
        timer=timer+1;
       end
    end
end
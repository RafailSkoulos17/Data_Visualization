
  var width = 1600;
  var height = 500;
  var padding = 60;
  var barWidth = 20;

  var margin = {top: 20, right: 20, bottom: 30, left: 50};

  var width = width - margin.left - margin.right,
      height = height - margin.top - margin.bottom;

  var totalWidth = width + margin.left + margin.right + padding;
  var totalheight = height + margin.top + margin.bottom + padding;
  
  d3.csv("team_salaries_2016-2017_extensively.csv", function(error, csv) {
    if (error) throw error;
    csv.sort(function(x,y) {return x.WinPer - y.WinPer});
    var data = [];
    var max = -Infinity;
    var min = +Infinity;
    csv.forEach(function(x) {
        var sal  = +x.Salary,
        t = x.Team,
        w = +x.WinPer,
        d = data[t];
        if (!d) {d = data[t] = [sal, w];}
        else d.push(sal)
        if (sal > max) max = sal;
        if (sal < min) min = sal;
    });
  // Sort group counts so quantile methods work
  for(var key in data) {
    var d = data[key];
    data[key] = d.sort(sortNumber);
  }
  var colorScale = d3.scaleOrdinal()
  .domain([1,30])
  .range(["#d684cd","#57c471","#b6439a","#96b23d","#746dd8","#d1972c","#4d328a",
  "#91bd68","#c771d2","#46b57c","#cd417f","#43c8ac","#d24b42","#5486e1",
  "#bda850","#4b5ba5","#cf733a","#66a1e5","#883016","#a787d4","#467528",
  "#6d2a70","#a17a37","#aa4f85","#ce7058","#822452","#d55a62","#e37ca7",
  "#8e2438","#d45e75"]);

  // Prepare the data for the box plots
  var boxPlotData = [];
  var i = 0;
  for (var [key, d] of Object.entries(data)) {
    d = d.map(Number);
    var record = {};
    var localMin = d3.min(d.slice(1));
    var localMax = d3.max(d.slice(1));
    record["rank"] = 30-i;
    record["key"] = key;
    record["winper"] = d.slice(0,1);
    record["counts"] = d.slice(1);
    record["quartile"] = boxQuartiles(d.slice(1));
    record["whiskers"] = [localMin, localMax];
    record["color"] = colorScale(i);
    i+= 1;
    boxPlotData.push(record);
  }
  // Compute an ordinal xScale for the keys in boxPlotData
  var xScale = d3.scalePoint()
    .domain(Object.keys(data))
    .rangeRound([padding, width - padding])
    .padding([1]);

  // Compute a global y scale based on the global counts
  var yScale = d3.scaleLinear()
    .domain([min-10, max*1.05])
    .range([height - padding, padding]);

  // Setup the svg and group we will draw the box plot in
  var svg = d3.select("#chart3").append("svg")
    .attr("width", totalWidth)
    .attr("height", totalheight)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // Setup the group the box plot elements will render in
  var g = svg.append("g")

    svg.append("text")
        .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
        .attr("transform", "translate("+ 1 +","+(height/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
        .attr('stroke-width', '5')
        .attr('font-size', 20)
        .text("Salary");


    svg.append("text")
    .attr("x", (width / 2))             
    .attr("y", 0 - ((margin.top / 2)-50))
    .attr("text-anchor", "middle")  
    .style("font-size", "26px") 
    .text("Salary distribution per Team (2016-2017)");

  // Draw the box plot vertical lines
  var verticalLines = g.selectAll(".verticalLines")
    .data(boxPlotData)
    .enter()
    .append("line")
    .attr("x1", function(datum) { return xScale(datum.key); })
    .attr("y1", function(datum) {
        var whisker = datum.whiskers[0];
        return yScale(whisker);
      }
    )
    .attr("x2", function(datum) { return xScale(datum.key); })
    .attr("y2", function(datum) {
        var whisker = datum.whiskers[1];
        return yScale(whisker);
      }
    )
    .attr("stroke", "#000")
    .attr("stroke-width", 1)
    .attr("fill", "none");

  // Draw the boxes of the box plot, filled in white and on top of vertical lines
  var rects = g.selectAll("rect")
    .data(boxPlotData)
    .enter()
    .append("rect")
    .attr("width", barWidth)
    .attr("height", function(datum) {
        var quartiles = datum.quartile;
        var height = yScale(quartiles[0]) - yScale(quartiles[2]);
        return height;
      }
    )
    .attr("x", function(datum) {
        return xScale(datum.key) - (barWidth/2) ;
      }
    )
    .attr("y", function(datum) {
        return yScale(datum.quartile[2]);
      }
    )
    .attr("fill", function(datum) {
      return datum.color;
      }
    )
    .attr("stroke", "#000")
    .attr("stroke-width", 1)
    .append("svg:title")
    .text(d => ` Team: ${d.key} 
Mean Salary: ${(d3.mean(d.counts)).toFixed(0)} 
Win Percentage: ${d.winper}
Rank: ${d.rank}`);
   

  // Now render all the horizontal lines at once - the whiskers and the median
  var horizontalLineConfigs = [
    // Top whisker
    {
      x1: function(datum) { return xScale(datum.key) - barWidth/2 },
      y1: function(datum) { return yScale(datum.whiskers[0]) },
      x2: function(datum) { return xScale(datum.key) + barWidth/2 },
      y2: function(datum) { return yScale(datum.whiskers[0]) }
    },
    // Median line
    {
      x1: function(datum) { return xScale(datum.key) - barWidth/2},
      y1: function(datum) { return yScale(datum.quartile[1]) },
      x2: function(datum) { return xScale(datum.key) + barWidth/2 },
      y2: function(datum) { return yScale(datum.quartile[1]) }
    },
    // Bottom whisker
    {
      x1: function(datum) { return xScale(datum.key) - barWidth/2 },
      y1: function(datum) { return yScale(datum.whiskers[1]) },
      x2: function(datum) { return xScale(datum.key) + barWidth/2 },
      y2: function(datum) { return yScale(datum.whiskers[1]) }
    }
  ];

  for(var i=0; i < horizontalLineConfigs.length; i++) {
    var lineConfig = horizontalLineConfigs[i];

    // Draw the whiskers at the min for this series
    var horizontalLine = g.selectAll(".whiskers")
      .data(boxPlotData)
      .enter()
      .append("line")
      .attr("x1", lineConfig.x1)
      .attr("y1", lineConfig.y1)
      .attr("x2", lineConfig.x2)
      .attr("y2", lineConfig.y2)
      .attr("stroke", "#000")
      .attr("stroke-width", 1)
      .attr("fill", "none");
  }

   //x-axis
   svg.append("g")
     .attr("transform", "translate(0," + (height - padding) + ")")
     .call(d3.axisBottom(xScale));

  // Add the Y Axis
  svg.append("g")
     .attr("transform", "translate("+ padding +",0)")
     .call(d3.axisLeft(yScale));

    });

	function boxQuartiles(d) {
  	return [
    	d3.quantile(d, .25),
    	d3.quantile(d, .5),
    	d3.quantile(d, .75)
  	];}

    
  // Perform a numeric sort on an array
  function sortNumber(a,b) {
    return a - b;
  }
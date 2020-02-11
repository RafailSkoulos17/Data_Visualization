d3.csv('team_stats_2007-2017.csv',function (data) {
    // CSV section
      // var body = d3.selectAll('body');
      data.Salary = +data.Salary;
      data.Efficiency = +data.Efficiency;
      data.WinPer = +data.WinPer;
      var active_link = "0";
      var active_team = "";
      var body = d3.select("#chart1");
      var formatPercent = d3.format('.1%');
      var formatThousand = d3.format(".0f");
      var selectData = [ { "text" : "Efficiency" },
                         { "text" : "Salary" }];
      var years = [ { "text" : "2017" },
                    { "text" : "2016" },
                    { "text" : "2015" },
                    { "text" : "2014" },
                    { "text" : "2013" },
                    { "text" : "2012" },
                    { "text" : "2011" },
                    { "text" : "2010" },
                    { "text" : "2009" },
                    { "text" : "2008" },
                    { "text" : "2007" }];

      // Select Y-axis Variable
      var span = body.append('span')
          .text('Select Y-Axis variable: ')

      var yInput = body.append('select')
          .attr('id','ySelect')
          .on('change',yChange)
        .selectAll('option')
          .data(selectData)
          .enter()
        .append('option')
          .attr('value', function (d) { return d.text })
          .text(function (d) { return d.text ;})
      body.append('br')

        
       // Select the year
      var span = body.append('span')
      .text('Select Year: ')

      var yearMenu = body.append("select")
          .attr('id','yearSelect')
          .on('change',yearChange)
        .selectAll("option")
          .data(years)
          .enter()
        .append("option")
          .attr("value", function(d){return d.text})
          .text(function(d){return d.text;})
    body.append('br')

    
      // Variables
      var body = d3.select('body');
      var margin = { top: 50, right: 50, bottom: 50, left: 50 };
      var padding = 80;
      var height = 500 - margin.top - margin.bottom;
      var width = 900 - margin.left - margin.right;
      var formatPercent = d3.format('.1%');
      
        var colorScale = d3.scaleOrdinal()
        .domain([1,30])
        .range(["#d684cd","#57c471","#b6439a","#96b23d","#746dd8","#d1972c","#4d328a",
        "#91bd68","#c771d2","#46b57c","#cd417f","#43c8ac","#d24b42","#5486e1",
        "#bda850","#4b5ba5","#cf733a","#66a1e5","#883016","#a787d4","#467528",
        "#6d2a70","#a17a37","#aa4f85","#ce7058","#822452","#d55a62","#e37ca7",
        "#8e2438","#d45e75"]);



      var xScale = d3.scaleLinear()
        .domain([
          d3.min(data,function (d) { return d['WinPer']*0.95 }),
          d3.max([0,d3.max(data,function (d) { return d['WinPer']*1.05 })])
          ])
        .range([padding,width - padding])
      var yScale = d3.scaleLinear()
        .domain([
          d3.min(data,function (d) { return d['Efficiency'] }),
          d3.max([0,d3.max(data,function (d) { return d['Efficiency'] *1.1})]) 
          ])
        .range([height - padding,padding])

      // SVG
      var svg = d3.select("#chart1").append('svg')
          .attr('height',height + margin.top + margin.bottom)
          .attr('width',width + margin.left + margin.right)
        .append('g')
          .attr('transform','translate(' + margin.left + ',' + margin.top + ')');
      // X-axis
      var xAxis = d3.axisBottom()
        .scale(xScale)
        .tickFormat(formatPercent);
        
      // Y-axis
      var yAxis = d3.axisLeft()
        .scale(yScale);
     

      
      // X-axis
      svg.append('g')
          .attr('class','axis')
          .attr('id','xAxis')
          .attr('transform', 'translate(0,' + (height - padding) + ')')
          .call(xAxis)
      
      svg.append("text")
          .attr('id','xAxisLabel')          
          .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
          .attr("transform", "translate("+ (width/2) +","+ height +")")  // centre below axis          .attr('stroke-width', '5')
          .attr('stroke-width', '5')
          .attr('font-size', 20)
          .text("Win Percentage");

      // Y-axis
      svg.append('g')
          .attr('class','axis')
          .attr('id','yAxis')
          .attr("transform", "translate("+ padding +",0)")
          .call(yAxis)
        

      svg.append("text")
          .attr('id', 'yAxisLabel')
          .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
          .attr("transform", "translate("+ (padding/20) +","+(height/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
          .attr('stroke-width', '5')
          .attr('font-size', 20)
          .text("Efficiency");
      
      svg.append("text")
      .attr("x", (width / 2))             
      .attr("y", 0 - ((margin.top / 2)-50))
      .attr("text-anchor", "middle")  
      .style("font-size", "26px") 
      .text("Efficiency and salary of teams over the years");
        
       // Create a dropdown
    
    var initializeGraph = function(year){

        // Filter the data to include only fruit of interest
        var selectYear = data.filter(function(d){
            return d.Year == year;
            })
            selectYear.sort(function(x,y) {return y.WinPer - x.WinPer});

        // Circles
        var circles = svg.selectAll('circle')
            .data(selectYear)
            .enter()
            .append('circle')
            .attr("class", "scatter")
            .attr('cx',function (d) { return xScale(d['WinPer']) })
            .attr('cy',function (d) { return yScale(d['Efficiency']) })
            .attr('r',10)
            .attr('stroke','black')
            .attr('stroke-width',1)
            .attr('opacity', 0.75)
            .attr('fill',function (d,i) { return colorScale(i) })

            .on('click', function (d){
                if (active_link == "0") {
                  active_link = "1"; 
                  active_team = d.Team; 
                  make_line_chart(d.Team);
                }
                else {
                  if (active_team == d.Team){
                    active_link = "0"; 
                    delete_line_chart();
                  }
                  else {
                    active_team = d.Team;
                    update_line_chart(d.Team);
                }
              }         
             } )
            .on('mouseover', function () {
                d3.select(this)
                .transition()
                .duration(500)
                .attr('r',20)
                .attr('stroke-width',3)
            })
            .on('mouseout', function () {
                d3.select(this)
                .transition()
                .duration(500)
                .attr('r',10)
                .attr('stroke-width',1)
            })
            .append('title') // Tooltip
            .text(function (d) { return 'Team: ' + d.Team +
                                '\nEfficiency: ' + (d['Efficiency']/1.0).toFixed(1) +
                                '\nWin Percentage: ' + formatPercent(d['WinPer']) +
                                '\nSalary: ' + (d['Salary']/1000000).toFixed(2) + " million dollars"});
        }

    initializeGraph("2017")
    console.log("graph initialized");
        // Update the data
 	var updateGraph = function(year){
    
    // find the previous choice of the user for y axis
    var yPrevChoice = d3.select("#ySelect").property("value");
    console.log(yPrevChoice);

    // Filter the data to include only year of interest
    var selectYear = data.filter(function(d){
            return d.Year == year;
            })
    selectYear.sort(function(x,y) {return y.WinPer - x.WinPer});
     // change the yScale
        yScale 
        .domain([
          d3.min([0,d3.min(selectYear,function (d) { return d[yPrevChoice] })]),
          d3.max([0,d3.max(selectYear,function (d) { return d[yPrevChoice]*1.1})])
          ])
      yAxis.scale(yScale) // change the yScale
      d3.select('#yAxis') // redraw the yAxis
        .transition().duration(1000)
        .call(yAxis)   

    xScale
        .domain([
          d3.min(selectYear,function (d) { return d['WinPer']*0.9 }),
          d3.max([0,d3.max(selectYear,function (d) { return d['WinPer'] *1.15})])
          ])

          xAxis.scale(xScale) // change the xScale
        d3.select('#xAxis') // redraw the xAxis
            .transition().duration(1000)
            .call(xAxis)   

           // Select all the lines and transition to new positions
          var selectYearGroups = svg.selectAll("circle")
            .data(selectYear);

          selectYearGroups.transition()
            .duration(1000)
            .delay(function (d,i) { return i*100})
              .attr('cy',function (d) { return yScale(d[yPrevChoice]) })
              .attr('cx',function (d) { return xScale(d['WinPer']) })
              .select("title")
              .text(function (d) { return 'Team: ' + d.Team +
                                '\nEfficiency: ' + (d['Efficiency']/1.0).toFixed(1) +
                                '\nWin Percentage: ' + formatPercent(d['WinPer']) +
                                '\nSalary: ' + (d['Salary']/1000000).toFixed(2) + " million dollars"})
    }

    function yearChange() {
       
        // get the selected year
        selectedYear = this.value;

       // update the graph for the selected year
       updateGraph(selectedYear);

}
      function yChange() {
        // var value = this.value // get the new y value
        var value = d3.select("#ySelect").property("value");

        console.log("SCATTER", value);
        yScale // change the yScale
          .domain([
            d3.min([0,d3.min(data,function (d) { return d[value] })]),
            d3.max([0,d3.max(data,function (d) { return d[value]*1.1})])
            ])
        yAxis.scale(yScale) // change the yScale
        d3.select('#yAxis') // redraw the yAxis
          .transition().duration(1000)
          .call(yAxis)
        d3.select('#yAxisLabel') // change the yAxisLabel
          .text(value)    
        d3.selectAll('circle.scatter') // move the circles
          .transition().duration(1000)
          .delay(function (d,i) { return i*100})
            .attr('cy',function (d) { return yScale(d[value]) })
      }


// -----------------------line-----------------------------------------

function delete_line_chart(){ d3.selectAll("#lineSVG").remove();
                              d3.selectAll("#lineSelect").remove();
                              d3.selectAll("#selText").remove();
                              d3.select("#chart2").selectAll("br").remove();
}


function update_line_chart(team) {

  var selectTeam = data.filter(function(d){
      return d.Team == team;
      })
  var value = d3.select("#lineSelect").property("value");

  var xScaleline = d3.scaleLinear()
      .domain([
          d3.min(data,function (d) { return +d['Year'] - 1 }),
          d3.max([0,d3.max(data,function (d) { return +d['Year'] + 1})])
          ])
      .range([padding,width - padding])
  
  
  var yScaleline = d3.scaleLinear()
      .domain([
          d3.min(data,function (d) { return +d[value] }),
          d3.max([0,d3.max(data,function (d) { return +d[value] })])
          ])
      .range([height - padding,padding]) 
  
  var radiusScaleline = d3.scaleLinear()
      .domain([d3.min(data, function(d){return +d.Salary;}), 
               d3.max([0, d3.max(data, function(d){return +d.Salary;})])])
      
      .range([3, 30]);
  
  // 7. d3's line generator
  var line = d3.line()
      .x(function(d) { return xScaleline(d.Year); }) // set the x values for the line generator
      .y(function(d) { return yScaleline(+d[value]); }) // set the y values for the line generator 
      .curve(d3.curveMonotoneX) // apply smoothing to the line
  
  my_svg = d3.selectAll("#lineSVG");
      my_svg.select("#line_yaxis").transition().duration(500).delay(function (d,i) { return i*5})
      .call(d3.axisLeft(yScaleline));
  
  d3.selectAll("path.line")
      .datum(selectTeam) // 10. Binds data to the line 
    .transition().duration(1000)
    .delay(function (d,i) { return i*100})
      .attr("class", "line") // Assign a class for styling 
      .attr("d", line); // 11. Calls the line generator 

// 12. Appends a circle for each datapoint 
d3.selectAll("circle.dot")
      .data(selectTeam)
      .transition()
      .duration(1000)
      .delay(function (d,i) { return i*5})// Uses the enter().append() method
          .attr("class", "dot") // Assign a class for styling
          .attr("cx", function(d) { return xScaleline(+d.Year) })
          .attr("cy", function(d) { return yScaleline(+d[value]) })
          .attr("r", 5)
          .select('title') // Tooltip
              .text(function (d) {return 'Team: ' + team +
                                  '\n Efficiency: ' + (d['Efficiency']/1.0).toFixed(1) +
                                  '\nWin Percentage: ' + formatPercent(d['WinPer']) +
                                  '\nSalary: ' + (+d['Salary']/1000000).toFixed(2) + " million dollars"});
}



function make_line_chart(team){

    var lineDataSelect = [ { "text" : "WinPer" },
                           { "text" : "Salary" },
                           { "text" : "Efficiency" }];

    var span = d3.select('#chart2').append('span')
          .attr("id", "selText")
          .text('Select Y-axis value: ')

    var yInput = d3.select('#chart2').append('select')
    .attr('id','lineSelect')
    .on('change',function() { lineChange(team); })
    .selectAll('option')
    .data(lineDataSelect)
    .enter()
    .append('option')
    .attr('value', function (d) { return d.text })
    .text(function (d) { return d.text;})
    d3.select('#chart2').append('br').append('br').append('br').append('br').append('br')

    var xScaleline = d3.scaleLinear()
        .domain([
            d3.min(data,function (d) { return +d['Year'] - 1 }),
            d3.max([0,d3.max(data,function (d) { return +d['Year'] + 1})])
            ])
        .range([padding,width - padding])
    
    
    var yScaleline = d3.scaleLinear()
        .domain([
            d3.min(data,function (d) { return d['WinPer'] })*0.9,
            d3.max([0,d3.max(data,function (d) { return d['WinPer'] })]) *1.1
            ])
        .range([height - padding,padding]) 
    
    var radiusScaleline = d3.scaleLinear()
        .domain([d3.min(data, function(d){return +d.Salary;}), 
                 d3.max([0, d3.max(data, function(d){return +d.Salary;})])])
        
        .range([3, 30]);
    
    
    // 7. d3's line generator
    var line = d3.line()
        .x(function(d) { return xScaleline(d.Year); }) // set the x values for the line generator
        .y(function(d) { return yScaleline(d.WinPer); }) // set the y values for the line generator 
        .curve(d3.curveMonotoneX) // apply smoothing to the line
    
    
    // 1. Add the SVG to the page and employ #2
    var line_svg = d3.select("#chart2").append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("id", "lineSVG")
        .attr("height", height + margin.top + margin.bottom)
      .append("g")
        .attr("transform", "translate(" + margin.left + "," + (margin.top) + ")");
    
    // 3. Call the x axis in a group tag
    line_svg.append("g")
        .attr("id", "line_xaxis")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + (height - padding) + ")")
        .call(d3.axisBottom(xScaleline).tickFormat(formatThousand)); // Create an axis component with d3.axisBottom
    
    // 4. Call the y axis in a group tag
    line_svg.append("g")
        .attr("id", "line_yaxis")
        .attr("class", "y axis")
        .attr("transform", "translate("+ padding +",0)")
        .call(d3.axisLeft(yScaleline)); // Create an axis component with d3.axisLeft
    
    line_svg.append("text")
        .attr('id', 'yAxisLabelLine')
        .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
        .attr("transform", "translate("+ (padding/20) +","+(height/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
        .attr('stroke-width', '5')
        .attr('font-size', 20)
        .text("Win Percentage");

  
    
    function initializeLine(team){
    
    // Filter the data to include only fruit of interest
    var selectTeam= data.filter(function(d){
                return d.Team == team;
                })
    

    line_svg.append("text")
    .attr("x", (width / 2))             
    .attr("y", 0 - ((margin.top / 2)-50))
    .attr("text-anchor", "middle")  
    .style("font-size", "26px") 
    .text("Statistics for " + team +" over the years");

    // 9. Append the path, bind the data, and call the line generator 
    line_svg.append("path")
        .datum(selectTeam) // 10. Binds data to the line 
        .attr("class", "line") // Assign a class for styling 
        .attr("d", line); // 11. Calls the line generator 
    
    // 12. Appends a circle for each datapoint 
    line_svg.selectAll("circles")
        .data(selectTeam)
      .enter().append("circle") // Uses the enter().append() method
        .attr("class", "dot") // Assign a class for styling
        .attr("cx", function(d) { return xScaleline(+d.Year) })
        .attr("cy", function(d) { return yScaleline(d.WinPer) })
        .attr("r", 5)
        .on('mouseover', function () {
                    d3.select(this)
                    .transition()
                    .duration(500)
                    .attr('r',15)
                    .attr('stroke-width',3)
                })
                .on('mouseout', function () {
                    d3.select(this)
                    .transition()
                    .duration(500)
                    .attr('r',5)
                    .attr('stroke-width',1)
                })
        .append('title') // Tooltip
            .text(function (d) { return 'Team: ' + team +
                                '\n Efficiency: ' + (d['Efficiency']/1.0).toFixed(1) +
                                '\nWin Percentage: ' + formatPercent(d['WinPer']) +
                                '\nSalary: ' + (d['Salary']/1000000).toFixed(2) + " million dollars"});
    
    };
    
    initializeLine(team);
    
      }

      function lineChange(team) {
      
      
      var value = d3.select("#lineSelect").property("value");
      var selectTeam = data.filter(function(d){
        return d.Team == team;
        })

      var line = d3.line()
          .x(function(d) { return xScaleline(d.Year); }) // set the x values for the line generator
          .y(function(d) { return yScaleline(d[value]); }) // set the y values for the line generator 
          .curve(d3.curveMonotoneX) // apply smoothing to the line

      var radiusScaleline = d3.scaleLinear()
      .domain([d3.min(data, function(d){return +d.Salary;}), 
               d3.max([0, d3.max(data, function(d){return +d.Salary;})])])
      .range([3, 30]);

      var xScaleline = d3.scaleLinear()
          .domain([
              d3.min(data,function (d) { return +d['Year'] - 1 }),
              d3.max([0,d3.max(data,function (d) { return +d['Year'] + 1})])
              ])
          .range([padding,width - padding])

      var yScaleline = d3.scaleLinear()
          .domain([
            d3.min([0,d3.min(data,function (d) { return d[value]*0.9 })]),
            d3.max([0,d3.max(data,function (d) { return d[value]*1.1})])
            ])
          .range([height - padding,padding]) 

      d3.select('#yAxisLabelLine') // change the yAxisLabel
      .text(value) 
      my_svg = d3.selectAll("#lineSVG");
      my_svg.select("#line_yaxis").transition().duration(500).delay(function (d,i) { return i*5})
      .call(d3.axisLeft(yScaleline));
  
      my_svg.selectAll("path.line")
          .datum(selectTeam) // 10. Binds data to the line 
          .transition().duration(1000)
          .delay(function (d,i) { return i*100})
          .attr("class", "line") // Assign a class for styling 
          .attr("d", line); // 11. Calls the line generator 

// 12. Appends a circle for each datapoint 
      my_svg.selectAll("circle")
          .data(selectTeam)
          .transition()
          .duration(1000)
          .delay(function (d,i) { return i*5})// Uses the enter().append() method
              .attr("class", "dot") // Assign a class for styling
              .attr("cx", function(d) { return xScaleline(+d.Year) })
              .attr("cy", function(d) { return yScaleline(d[value]) })
              .attr("r", 5)
              .select('title') // Tooltip
                  .text(function (d) {return 'Team: ' + team +
                                      '\n Efficiency: ' + (d['Efficiency']/1.0).toFixed(1) +
                                      '\nWin Percentage: ' + formatPercent(d['WinPer']) +
                                      '\nSalary: ' + (d['Salary']/1000000).toFixed(2) + " million dollars"});
        }
      
    })
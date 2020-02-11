/*
              START
    >>>>TREE SELECTION TOOL<<<<
*/
d3.csv("player_stats_2007-2017.csv", function(data){

data = data.filter(function(x){return x.Year == "2017"}) //Salary is higher than 5 millions
data = data.filter(function(x){return x.Tm !== "TOT"}) //Salary is higher than 5 millions
data = data.filter(function(x){return x.Salary>5000000}) //Salary is higher than 5 millions
data = data.filter(function(x){return x.EffPerSalary*10000000>0}) //Salary is higher than 5 millions
var treeData= { "key": "NBA", "values": 
        d3.nest()
        .key(function (d) { return d.Tm; })
        .entries(data)
      };

//console.log(treeData)
var treeData = { "name": "NBA", "children":
        treeData.values.map( function(major){

            return { "name": major.key, "children": 
              major.values.map( function(region){


                      return { "name": region.Player };                  
                 
              }) 
            }; 

         }) 
      };

//console.log(treeData)

var margin = {top: 20, right: 20, bottom: 30, left: 20},
    width = 1800 - margin.left - margin.right,
    height = 300 - margin.top - margin.bottom;

var svg = d3.select("#chart0")
        .append("svg")
        .attr("id","the_SVG_ID2")
        .attr("width", width + margin.right + margin.left)
        .attr("height", height + margin.top + margin.bottom)

var g3 = svg.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

// Setup SVG Element - End

var i = 0,
    duration = 750,
    root;

// Setup tree

var treemap = d3.tree()
    //.size([width, height])
    .separation(function(a, b) { return ((a.parent == root) && (b.parent == root)) ? 1 : 2.5; })
    .size([width, width - 160]);

// Get the root

root = d3.hierarchy(treeData, function(d) { return d.children; });

root.x0 = 0;
root.y0 = width / 3;

// Collapse all children, except root's
//root.children.forEach(collapse);
root.children.forEach(collapse);
// root.children = null;

// Let's draw the tree
draw(root);

// console.log(root);

function draw(source) {

  // Get the treemap, so that we can get nodes and links
  var treeData = treemap(root);

  // Get nodes and links
  var nodes = treeData.descendants()
      links = treeData.descendants().slice(1);

  // Adjust the position of y of each node. Comment out just this line and see how it's different  
  nodes.forEach(function(d){ d.y = d.depth * 70});

  // Add unique id for each node, else it won't work
  var node = g3.selectAll('g.node')
      .data(nodes, function(d) {return d.id || (d.id = ++i);   });
      


  // Let's append all enter nodes
  var nodeEnter = node
      .enter()
      .append('g')
      .attr('class', 'node')
      .attr("transform", function(d) {
        return "translate(" + source.x0 + "," + source.y0 + ")";
      })
      .on('click',click);

  // Add circle for each enter node, but keep the radius 0

  nodeEnter.append('circle')
      .attr('class', 'node')
      .attr('r', 1e-6)
      .style("fill", function(d) {
          return d._children ? "lightsteelblue" : "#fff";
      });

  // Add text

  nodeEnter.append('text')
    .attr("dy", ".35em")
    //.attr('x', 20)
    .attr('y', 20)
    .attr("x", function(d) {
        return d.children || d._children ? -13 : 13;
    })
    .attr("text-anchor", function(d) {
        return d.children || d._children ? "start" : "end";
    })

    .text(function(d) { return d.data.name; });


  var nodeUpdate = nodeEnter.merge(node);

  // Do transition of node to appropriate position
  nodeUpdate.transition()
    .duration(duration)
    .attr("transform", function(d) { 
        return "translate(" + d.x + "," + d.y + ")";
     });


  // Let's update the radius now, which was previously zero.

  nodeUpdate.select('circle.node')
    .attr('r', 15 )
    .style("fill", function(d) {
        return d._children ? "lightsteelblue" : "#fff";
    })
    .attr('cursor', 'pointer');

  // Let's work on exiting nodes
  
  // Remove the node
  
  var nodeExit = node.exit().transition()
      .duration(duration)
      .attr("transform", function(d) {
          return "translate(" + source.x + "," + source.y + ")";
      })
      .remove();

  // On exit reduce the node circles size to 0
  nodeExit.select('circle')
    .attr('r', 1e-6);

  // On exit reduce the opacity of text labels
  nodeExit.select('text')
    .style('fill-opacity', 1e-6);

  
  // Let's draw links

  var link = g3.selectAll('path.link')
      .data(links, function(d) { return d.id; });
  
  // Work on enter links, draw straight lines

  var linkEnter = link.enter().insert('path', "g")
      .attr("class", "link")
      .attr('d', function(d){
        var o = {x: source.x0, y: source.y0}
        return diagonal(o, o)
      });

  // UPDATE
  var linkUpdate = linkEnter.merge(link);

  // Transition back to the parent element position, now draw a link from node to it's parent
  linkUpdate.transition()
      .duration(duration)
      .attr('d', function(d){ return diagonal(d, d.parent) });

  // Remove any exiting links
  var linkExit = link.exit().transition()
      .duration(duration)
      .attr('d', function(d) {
        var o = {x: source.x, y: source.y}
        return diagonal(o, o)
      })
      .remove();

  // Store the old positions for transition.
  nodes.forEach(function(d){
    d.x0 = d.x;
    d.y0 = d.y;
  });

}
  
function diagonal(s, d) {
  

  
   var path = `M ${s.x} ${s.y}
           C ${(s.x + d.x) / 2} ${s.y},
             ${(s.x + d.x) / 2} ${d.y},
             ${d.x} ${d.y}`

  return path
}

function collapse(d) {
  if(d.children) {
    d._children = d.children
    d._children.forEach(collapse)
    d.children = null
  }
}

function click(d)
{
  var click = d.data.name
  console.log(click)

  if(click.length < 4 ) {
    console.log("Your name must contain minimum 4 characters");
    // return false;
    change(click)
  }
  else{stacked_bar(click)}


  if (d.children) {
      d._children = d.children;
      d.children = null;
    } else {
      d.children = d._children;
      d._children = null;
    }
  // If d has a parent, collapse other children of that parent
  if (d.parent) {
    d.parent.children.forEach(function(element) {
      if (d !== element) {
        collapse(element);
      }
    });
  }

  draw(d);
  //change();
}
 //});
/* 
              END
    >>>>TREE SELECTION TOOL<<<<
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
                 START
    >>>BAR CHART ADN STACKED BAR CHART<<<
*/

/*
 Bar chart for the players of each team, for each player we visualize the efficiency per
 salary to see who are the most value for money players, we visualize players with positive
 efficiency and players with salary that is higher than 5 million.
 */

margin = ({top: 20, right: 20, bottom: 20, left: 60})
height = 400
width = 800
var padding = 70;
var active_link = "0"; //to control legend selections and hover
var legendClicked; //to control legend selections
var legendClassArray = []; //store legend classes to select bars in plotSingle()
var y_orig; //to store original y-posn
//function chart(){
//d3.csv("bar_chart_2018.csv", function(data){


/*
We make an array of arrays the keys are the teams as we can see from the console.log()
*/
var dataNested = d3.nest()
  .key(function (d) { return d.Tm; })
  .entries(data)


function change(a) {
    d3.select("#identity").remove();
    restorePlot();

    value = a

    
    dataFiltered = dataNested.filter(function (d) { return d.key === value})
    console.log(value)
    dataFiltered = dataFiltered[0].values


    dataFiltered = dataFiltered.filter(function(x){return x.Salary>5000000}) //Salary is higher than 5 millions
    dataFiltered = dataFiltered.filter(function(d){return d.EffPerSalary>0}) //We dont visualize players that have negative efficiency

    var dataFiltered = dataFiltered.map(function(d,i) {

      return {
        Player: d.Player,
        TM: d.Tm,
        Salary: d.Salary,
        EffPerSalary: d.EffPerSalary*10000000
      };
    });


    var efficiency_avg = d3.mean(dataFiltered, function(d) { return d.EffPerSalary; }); // the average efficinecy


    x = d3.scaleBand()
    .domain(dataFiltered.map(function(d){return d.Player}))
    .range([margin.left, width - margin.right])
    .padding(0.1)


    y = d3.scaleLinear()
    .domain([0, d3.max(dataFiltered, d => d.EffPerSalary*1)]).nice() //Multiply it by 1 (unjustified)???
    .range([height - margin.bottom, margin.top])


    xAxis = g => g
    .attr("transform", `translate(0,${height - margin.bottom})`)
    .call(d3.axisBottom(x)
    .tickSizeOuter(0))
    

    yAxis = g => g
      .attr("transform", `translate(${margin.left},0)`)
      .call(d3.axisLeft(y))
      .call(g => g.select(".domain").remove())

      var svg = d3.select("#chart0")
            .append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
            .attr("id", "identity" )
      
          svg.append("text")
            .attr('id', 'yAxisLabel')
            .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
            .attr("transform", "translate("+ ((padding+400)/20) +","+(height/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
            .attr('stroke-width', '5')
            .attr('font-size', 20)
            .text("Efficiency per Salary");
            


          svg.append("text")
            .attr("x", (width / 2))             
            .attr("y", 0 - ((margin.top)-50))
            .attr("text-anchor", "middle")  
            .style("font-size", "26px") 
            .text("Efficiency per salary for the players of team:"+ value);


      svg.append("g")
      .attr("fill", "#9a5ea1")
      .style("opacity", "0.5")
      .selectAll("rect").data(dataFiltered).enter().append("rect")

      .attr("x", function(d) {return x(d.Player);}) 
      .attr("y", d => y(d.EffPerSalary))
      .attr("height", d => y(0) - y(d.EffPerSalary))  
      .attr("width", x.bandwidth())
      .on('mouseover', function () {
        d3.select(this)
        .style('fill',"#98823c")
        .style('stroke','black')
        .style('stroke-width',3)
    })
    .on('mouseout', function () {
      d3.select(this)
      .style("fill", "#9a5ea1")
      if (active_link == "0"){
        d3.select(this)
          .style('stroke','none')
      }
      
  })
      .on("click",function(d){        

      if (active_link === "0") //nothing  , turn on this selection
      { 
          d3.select(this)           
            .style("stroke", "black")
            .style("stroke-width", 2);

          active_link = this.id.split("id").pop();
          
          var value2 = d.Player
          stacked_bar(value2);  //make the stacked_bar_graph

          //gray out the others
          for (i = 0; i < legendClassArray.length; i++)
          {
            if (legendClassArray[i] != active_link) 
            {
              d3.select("#id" + legendClassArray[i])
                .style("opacity", 0.5);
            }
          }
          
          } 
      else { //deactivate
        if (active_link === this.id.split("id").pop()) //active square selected; turn it OFF
        {
          d3.select(this)           
            .style("stroke", "none");

          active_link = "0"; //reset

          //restore remaining boxes to normal opacity
          for (i = 0; i < legendClassArray.length; i++) 
          {              
              d3.select("#id" + legendClassArray[i])
                .style("opacity", 1);
          }

              //restore plot to original
              
              restorePlot();

        }
            
    
            } //end active_link check
            
                              
                                    
          });
          
          svg.append("line")
          .attr("x1", margin.left)
          .attr("y1", y(efficiency_avg))
          .attr("x2", width)
          .attr("y2", y(efficiency_avg))
          .attr("stroke-width", 2)
          .attr("stroke", "#ee8866");

          svg.append("g")
          .call(xAxis);
      
          svg.append("g")
            .call(yAxis);
  }
  });
//}
/*
Function that should restore the Plot but now just shows a second graph under the first graph
*/

function restorePlot(d) {
  d3.selectAll("#the_SVG_ID").remove();
  
    }

//chart();



/*
>>>>Function that shows the stacked_bar graph<<<<
*/
//

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function stacked_bar(value2){
  restorePlot();
  
  margin = ({top: 20, right: 20, bottom: 30, left: 60})
  height = 400
  width = 700
  padding = 70
    
    var svg2 = d3.select("#chart0")
    .append("svg")
    .attr("id","the_SVG_ID")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)

    svg2.append("text")
    .attr('id', 'yAxisLabel')
    .attr("text-anchor", "middle")  // this makes it easy to centre the text as the transform is applied to the anchor
    .attr("transform", "translate("+ ((padding+400)/20) +","+(height/2)+")rotate(-90)")  // text is drawn off the screen top left, move down and out and rotate
    .attr('stroke-width', '5')
    .attr('font-size', 20)
    .text("Total #accomplishments");   


    svg2.append("text")
    .attr("x", (width / 1.5))             
    .attr("y", 0 - ((margin.top)-50))
    .attr("text-anchor", "middle")  
    .style("font-size", "26px") 
    .text("Accomplishments of " + value2 + " over the years");
    // .style("opacity", "0.5");

    g2 = svg2.append("g").attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    
    // set x scale
    var x = d3.scaleBand()
        .rangeRound([0, width])
        .paddingInner(0.05)
        .align(0.1);
    
    // set y scale
    var y = d3.scaleLinear()
        .rangeRound([height, 0]);
    
    // set the colors
    var z = d3.scaleOrdinal()
    
    .range(['#77aadd', '#88ddff', '#44bb99', '#aaaa00', '#ee8866','#ffaabb']);
    



    d3.csv("player_stats_2007-2017.csv", function(data_bar) {

      var dataNested_Player = d3.nest()
      .key(function (d) { return d.Player; })
      .entries(data_bar)

    var data_bar = dataNested_Player.filter(function (d) { return d.key === value2 })
    console.log(value2)

    var data_LeBron = data_bar[0].values

    function round(num) {
      return +(Math.round(num + "e+2")  + "e-2");
    }
  
    var data_LeBron = data_LeBron.map(function(d) {

    

      return { 
        Team: d.Tm,
        Date: d.Year,
        Player: d.Player,
        Three_points: round(d["3P"]/d["G"]),
        Two_points: round(d["2P"]/d["G"]),
        Free_throws: round(d.FT/d["G"]),
        
        Offensive_rebounds: round(d.ORB/d["G"]),
        Deffensive_rebounds: round(d.DRB/d["G"]),
        
        Assists: round(d.AST/d["G"])
      };
    });
    
    /*
    array with Points, Rebounds, Assists
    */

   const test = data_LeBron.map((d, index, array) =>{return d.Date;
  });

  

   var dataNested_date = d3.nest()
   .key(function (d) { return d.Date; })
   .entries(data_LeBron)

  var data_test = [];
    for (var i = 0; i<dataNested_date.length; i++)
    {
     data_test[i] = dataNested_date[i].values[0];
    }




    var data_LeBron_total = data_test.map(function(d, i) {
    


       return { 
         Three_points: d.Three_points,
         Two_points: d.Two_points,
         Free_throws: d.Free_throws,
         Offensive_rebounds: d.Offensive_rebounds,
         Deffensive_rebounds: d.Deffensive_rebounds,

         Assists: d.Assists
       };
     });

    /*
    array with the sum of Points, Rebounds, Assists
    */
     
     var data_LeBron_total_number = data_LeBron_total.map(function(d) {
       
       function sum( obj ) {
        var sum = 0;
        for( var el in obj ) {
          if( obj.hasOwnProperty( el ) ) {
            sum += parseFloat( obj[el] );
          }
        }
        //console.log(sum)
        return sum;
      }
      //console.log(d)
      to = sum(d)
      //console.log(to)
        return {total: to };
    });




    //>>>>>>Add each stats for the player that changed teams

    keys = d3.keys(data_test[0])
    keys = keys.slice(-6);

  x.domain(data_LeBron.map(function(d) { return d.Date; }));
  y.domain([0, d3.max(data_LeBron_total_number, function(d) {
     return d.total; })]).nice();
  z.domain(keys);

  g2.append("g")
    .selectAll("g")
    .data(d3.stack().keys(keys)(data_test))
    .enter().append("g")
      .attr("fill", function(d) {
        return z(d.key); })
    .style("opacity", "0.5")
    .selectAll("rect")
    .data(function(d) { return d; })
    .enter().append("rect")
      .attr("x", function(d) {
         return x(d.data.Date); })
      .attr("y", function(d) { 
        return y(d[1]); })
      .attr("height", function(d) { return y(d[0]) - y(d[1]); })  
      .attr("width", x.bandwidth())
    .on("mouseover", function() { tooltip.style("display", null); })
    .on("mouseout", function() { tooltip.style("display", "none"); })
    .on("mousemove", function(d) {

      var xPosition = d3.mouse(this)[0] - 5;
      var yPosition = d3.mouse(this)[1] - 5;
      tooltip.attr("transform", "translate(" + xPosition + "," + yPosition + ")");
      tooltip.select("text").text(d[1]-d[0]);
    });
  g2.append("g")
      .attr("class", "axis")
      .attr("transform", "translate(0," + height + ")")
      .call(d3.axisBottom(x));

  g2.append("g")
      .attr("class", "axis")
      .call(d3.axisLeft(y).ticks(null, "s"))
    .append("text")
      .attr("x", 2)
      .attr("y", y(y.ticks().pop()) + 0.5)
      .attr("dy", "0.32em")
      .attr("fill", "#000")
      .attr("font-weight", "bold")
      .attr("text-anchor", "start");

  var legend = g2.append("g")
      .attr("font-family", "sans-serif")
      .style("opacity", "0.7")
      .attr("font-size", 10)
      .attr("text-anchor", "start")
    .selectAll("g")
    .data(keys.slice().reverse())
    .enter().append("g")
      .attr("transform", function(d, i) { return "translate(-675," + i * 20 + ")"; });

  legend.append("rect")
      .attr("x", width - 19)
      .attr("width", 19)
      .attr("height", 19)
      .attr("fill", z);

  legend.append("text")
      .attr("x", width + 4)
      .attr("y", 9.5)
      .attr("dy", "0.32em")
      .text(function(d) { return d; });
});

  // Prep the tooltip bits, initial display is hidden
  var tooltip = svg2.append("g")
    .attr("class", "tooltip")
    .style("display", "none");
      
  tooltip.append("rect")
    .attr("width", 60)
    .attr("height", 20)
    .attr("fill", "white")
    .style("opacity", 0.5);

  tooltip.append("text")
    .attr("x", 30)
    .attr("dy", "1.2em")
    .style("text-anchor", "middle")
    .attr("font-size", "12px")
    .attr("font-weight", "bold");

    }

    /*
                 END
    >>>BAR CHART ADN STACKED BAR CHART<<<
*/
package main

import (
	"math/rand"
  "fmt"
  "time"
  "strconv"
  "image"
  "image/png"
  "os"
  "bufio"
  "sort"
  "math"
  "io/ioutil"
  "github.com/wcharczuk/go-chart"
	"github.com/faiface/pixel"
	"github.com/faiface/pixel/imdraw"
	"github.com/faiface/pixel/pixelgl"
	"golang.org/x/image/colornames"
  "github.com/fatih/color"
  "github.com/disintegration/imaging"
  "github.com/manifoldco/promptui"
  //"github.com/k0kubun/pp"
)



const (
	maxpopulation  = 1000
  maxpolygons = 1000
	maxvertices = 50
)

var population int = 100
var polygons int = 100
var vertices int = 3
var tournamentsize int = 2
var elite int = 1
var cutoff int = 50
var genlog int = 10
var xsize int = 200
var ysize int = 246
var mutationrate float64 = 2.0

type Organism struct {
	DNA
	Fitness float64
}

type DNA struct {
  Genes [maxpolygons]Gene
}

type Gene struct{
  r float64
  g float64
  b float64
  a float64
  Points [maxvertices]Point
}

type Point struct{
  x int
  y int
}

var Population [maxpopulation]Organism



var Generation int = 0
var Fitnesses []float64 = make( []float64,0 )
var refimage []int
var currentGenePool []Organism = make( []Organism,0 )
var logpath string
var loggen int = 0


func createInitialPopulation() {

  for i := 0; i < population; i++ {

    var genes [maxpolygons]Gene

    for j := 0; j < polygons; j++ {

      var points [maxvertices]Point

      for k :=0; k < vertices; k++ {
        point := Point{
          x: rand.Intn(xsize),
          y: rand.Intn(ysize),
        }
        points[k] = point
      }

      gene := Gene{
        r: rand.Float64(),
        g: rand.Float64(),
        b: rand.Float64(),
        a: rand.Float64(),
        Points: points,
      }
      genes[j] = gene

    }

    dna := DNA{
      Genes: genes,
    }

    organism := Organism{
      DNA: dna,
      Fitness: 0,
    }

    Population[i] = organism
  }

}

func calculateFitness(){

  for i := 0; i < population; i++ {

    //Virtual canvas
    imd := *imdraw.New(nil)
    drawOrganism(Population[i], &imd)


    a := pixel.Vec {
      X: 0.0,
      Y: 0.0,
    }
    b := pixel.Vec {
      X: float64(xsize),
      Y: float64(ysize),
    }

    bounds := pixel.Rect {
      Min: a,
      Max: b,
    }

    canvas := pixelgl.NewCanvas(bounds)
    imd.Draw(canvas)

    pixels := canvas.Pixels()

    pixelsint := make([]int, len(pixels))

    for k:=0; k < len(pixels); k++{
      pixelsint[k] = int(pixels[k])
    }



    fitness := 0.0


    for j := 0; j < xsize * ysize * 4; j += 4 {

      deltaRed := pixelsint[j] - refimage[j]
      deltaGreen := pixelsint[j + 1] - refimage[j + 1]
      deltaBlue := pixelsint[j + 2] - refimage[j + 2]

      //Color-Distance in 3D-Space
      pixelFitness := float64(deltaRed * deltaRed + deltaGreen * deltaGreen + deltaBlue * deltaBlue)
      fitness -= pixelFitness

    }

    Population[i].Fitness = -1 * math.Sqrt(-1 * (fitness / 10000000));

  }

}

func bestFitness() (int, float64){

  sortPopulation()
  return 0, Population[0].Fitness

}

func averageFitness() (float64){

  var totalFitness float64

  for i := 0; i < population; i++ {

    totalFitness += Population[i].Fitness

  }

  totalFitness = totalFitness / float64(population)
  return totalFitness

}

/*func createGenePool(maxFitness float64){

  var newGenePool []Organism = make( []Organism,0 )

  for i := 0; i < population; i++ {

		num := int((Population[i].Fitness / maxFitness ) * 100)
		for n := 0; n < num; n++ {
			newGenePool = append(newGenePool, Population[i])
		}

	}

  currentGenePool = newGenePool


}*/

func TournamentSelection() Organism{

  var tournament []Organism = make( []Organism,0 )

  for n:=0; n < tournamentsize; n++{

    i := rand.Intn(population - cutoff)
    tournament = append(tournament, Population[i])

  }

  var bestint int = 0
  var bestfit float64 = tournament[0].Fitness

  for n:=0; n < tournamentsize; n++{

    if tournament[n].Fitness > bestfit {
      bestint = n
      bestfit = tournament[n].Fitness
    }

  }

  return tournament[bestint]

}

func sortPopulation () { //sort by highest fitness

  sorted := make([]Organism, population)

  for i :=0;i < population; i++{

    sorted[i] = Population[i]

  }

  sort.Slice(sorted, func(i, j int) bool {
  return sorted[i].Fitness > sorted[j].Fitness
  })

  for j :=0; j < population; j++{

    Population[j] = sorted[j]

  }

}

func naturalSelection() {
	var nextPopulation [maxpopulation]Organism

  sortPopulation()

	for i := 0; i < population; i++ {

    if i < elite {
      nextPopulation[i] = Population[i]

    }else {

      a := TournamentSelection()
      b := TournamentSelection()


		  child := Crossover(a, b)
		  child.Mutate()


		  nextPopulation[i] = child
    }

	}

  Population = nextPopulation
}

//1 point crossover
func Crossover(p1 Organism, p2 Organism) Organism {

  var genes [maxpolygons]Gene

  dna := DNA{
    Genes: genes,
  }

	child := Organism{
		DNA:    dna,
		Fitness: 0,
	}

	mid := rand.Intn(polygons)

	for i := 0; i < polygons; i++ {

		if i > mid {
			child.DNA.Genes[i] = p1.DNA.Genes[i]
		} else {
			child.DNA.Genes[i] = p2.DNA.Genes[i]
		}

	}

	return child
}


func (o *Organism) Mutate() {

	for i := 0; i < polygons; i++ {

		if rand.Float64() * 100.0 < mutationrate {

      var index int = rand.Intn(6) 

      if index == 0 {
        o.DNA.Genes[i].r = rand.Float64()
      }
      if index == 1 {
        o.DNA.Genes[i].g = rand.Float64()
      }
      if index == 2 {
        o.DNA.Genes[i].b = rand.Float64()
      }
      if index == 3 {
        o.DNA.Genes[i].a = rand.Float64()
      }
      if index == 4 {

        var gene Gene
        var j int = rand.Intn(polygons)

        gene = o.DNA.Genes[j]
        o.DNA.Genes[j] = o.DNA.Genes[i]
        o.DNA.Genes[i] = gene

      }
      if index > 4 {

        index = rand.Intn(vertices)

        point := Point{
          x: rand.Intn(xsize),
          y: rand.Intn(ysize),
        }
        o.DNA.Genes[i].Points[index] = point

      }

		}
	}

}


func drawOrganism( o Organism, imd *imdraw.IMDraw){

  //background
  imd.Color = pixel.RGB(1, 1, 1)
  imd.Push(pixel.V(0, 0))
  imd.Push(pixel.V(0, float64(ysize)))
  imd.Push(pixel.V(float64(xsize), float64(ysize)))
  imd.Push(pixel.V(float64(xsize), 0))
  imd.Polygon(0)
  ////

  //vertices
  for i:=0; i < polygons; i++{

    gene := o.DNA.Genes[i]

    color := pixel.RGB(gene.r, gene.g, gene.b).Mul(pixel.Alpha(gene.a))

    imd.Color = color

    for j:=0; j < vertices; j++{

      point := gene.Points[j]
      imd.Push(pixel.V(float64(point.x), float64(point.y)))

    }
    imd.Polygon(0)

  }

}


func run() {


  loadReferenceImage()
  printConfig()

  //adjust window size
  var ys int
  if (ysize < 400){
    ys = 400
  }else{
    ys = ysize
  }
  ///

  //setup window
	cfg := pixelgl.WindowConfig{
		Title:  "Genetic Algorithm",
		Bounds: pixel.R(0, 0, float64(xsize + 10 + 1024), float64(ys)),
		//VSync:  true,
	}
	win, err := pixelgl.NewWindow(cfg)

	if err != nil {
		panic(err)
	}
  //

  startLog()


  //initialize algorithm
  createInitialPopulation()
  calculateFitness()
  //

  //display best from first generation
  bestint, bestfit := bestFitness()
  imd := *imdraw.New(nil)
  drawOrganism(Population[bestint], &imd)
  Fitnesses = append(Fitnesses, bestfit)
  //

  //display info
  printGenerationInfo(bestint)

  //update window
  win.Clear(colornames.White)
  imd.Draw(win)
  win.Update()



	for !win.Closed() { // loop for each generation


    calculateFitness()
    bestint, bestfit := bestFitness()
    Fitnesses = append(Fitnesses, bestfit)


    naturalSelection()

    calculateFitness()
    bestint, bestfit = bestFitness()

    Generation += 1


    imd := *imdraw.New(nil)
    drawOrganism(Population[bestint], &imd)

    printGenerationInfo(bestint)

    //Draw chart
    pic, err := loadPicture("chart.png")
	  if err != nil {
      err = nil
	  }else{
	  sprite := pixel.NewSprite(pic, pic.Bounds())
    win.Clear(colornames.White)
    sprite.Draw(win, pixel.IM.Moved( pixel.V(float64(xsize + 10 + 512), (win.Bounds().Max.Y - win.Bounds().Min.Y) / 2)  ))
    }
    /////

    //update window
		imd.Draw(win)
		win.Update()

	}

}










//Helper functions


func startLog(){




  input, err := ioutil.ReadFile("ref.png")
       if err != nil {
               fmt.Println(err)
               return
       }

       err = ioutil.WriteFile(logpath + "ref.png", input, 0644)
       if err != nil {
               fmt.Println("Error creating", logpath + "ref.png")
               fmt.Println(err)
               return
       }




}

func logResult(bestint int) {

  //Virtual canvas
  imd := *imdraw.New(nil)
  drawOrganism(Population[bestint], &imd)


  a := pixel.Vec {
    X: 0.0,
    Y: 0.0,
  }
  b := pixel.Vec {
    X: float64(xsize),
    Y: float64(ysize),
  }

  bounds := pixel.Rect {
    Min: a,
    Max: b,
  }

  canvas := pixelgl.NewCanvas(bounds)
  imd.Draw(canvas)

  pixels := canvas.Pixels()



  imgf := image.NewRGBA(image.Rect(0, 0, xsize, ysize))
  imgf.Pix = pixels

  img := imaging.FlipV(imgf)

  savpath := logpath + "i" + string(strconv.Itoa(Generation)) + ".png"
  fmt.Println(savpath)

  f, _ := os.OpenFile(savpath, os.O_WRONLY|os.O_CREATE, 0600)
    defer f.Close()
  png.Encode(f, img)

  input, err := ioutil.ReadFile("chart.png")
       if err != nil {
               fmt.Println(err)
               return
       }

       err = ioutil.WriteFile(logpath + "p" + string(strconv.Itoa(Generation)) + ".png", input, 0644)
       if err != nil {
               fmt.Println(err)
               return
       }



}

func printGenerationInfo(bestint int) {

  cyan := color.New(color.FgCyan).PrintlnFunc()
  fmt.Printf("Current generation: ")
  cyan(strconv.Itoa(Generation))
  fmt.Printf("Best Fitness: ")
  cyan(Population[bestint].Fitness)
  fmt.Printf("Average Fitness: ")
  cyan(averageFitness())
  fmt.Println("")
  drawChart()

  if int(float64(Generation / genlog)) > loggen || Generation == 0  {
    loggen = int(float64(Generation / genlog))
    logResult(bestint)

  }


}

func printConfig(){

  cyan := color.New(color.FgCyan).PrintlnFunc()
  fmt.Printf("Initial population: ")
  cyan(strconv.Itoa(population))
  fmt.Printf("Genes per organism: ")
  cyan(strconv.Itoa(polygons))
  fmt.Printf("Vertices per gene: ")
  cyan(strconv.Itoa(vertices))
  fmt.Printf("Tournament size: ")
  cyan(tournamentsize)
  fmt.Printf("Elite: ")
  cyan(elite)
  fmt.Printf("Cutoff: ")
  cyan(cutoff)
  fmt.Printf("Mutationrate: ")
  cyan(mutationrate)
  fmt.Printf("Log every x generation: ")
  cyan(genlog)
  fmt.Println(">")
  fmt.Println(">")
  fmt.Println(">")
  fmt.Println(">")

  prompt := promptui.Prompt{
		Label:     "Use this config",
		IsConfirm: true,
	}

	result, _ := prompt.Run()

  if (result == "y"){

  }else{



  	prompt2 := promptui.Prompt{
  		Label:    "Population size",
  	}

  	result2, _ := prompt2.Run()

    population = stringint(result2)

    prompt2 = promptui.Prompt{
  		Label:    "Genes per Organism",
  	}

  	result2, _ = prompt2.Run()

    polygons = stringint(result2)

    prompt2 = promptui.Prompt{
  		Label:    "Vertices per Gene",
  	}

  	result2, _ = prompt2.Run()

    vertices = stringint(result2)

    prompt2 = promptui.Prompt{
  		Label:    "Tournament size",
  	}

  	result2, _ = prompt2.Run()

    tournamentsize = stringint(result2)

    prompt2 = promptui.Prompt{
  		Label:    "Amount of elite Organisms",
  	}

  	result2, _ = prompt2.Run()

    elite = stringint(result2)

    prompt2 = promptui.Prompt{
  		Label:    "Cutoff amount",
  	}

  	result2, _ = prompt2.Run()

    cutoff = stringint(result2)

    prompt2 = promptui.Prompt{
  		Label:    "Mutationrate (%) per organism",
  	}

  	result2, _ = prompt2.Run()

    mutationrate = float64(stringint(result2))

    prompt2 = promptui.Prompt{
  		Label:    "Log every x amount of generations",
  	}

  	result2, _ = prompt2.Run()

    genlog = stringint(result2)



    cyan := color.New(color.FgCyan).PrintlnFunc()
    fmt.Println(">")
    fmt.Println(">")
    fmt.Println(">")
    fmt.Println(">")
    fmt.Printf("Initial population: ")
    cyan(strconv.Itoa(population))
    fmt.Printf("Genes per organism: ")
    cyan(strconv.Itoa(polygons))
    fmt.Printf("Vertices per gene: ")
    cyan(strconv.Itoa(vertices))
    fmt.Printf("Tournament size: ")
    cyan(tournamentsize)
    fmt.Printf("Elite: ")
    cyan(elite)
    fmt.Printf("Cutoff: ")
    cyan(cutoff)
    fmt.Printf("Mutationrate: ")
    cyan(mutationrate)
    fmt.Printf("Log every x generation: ")
    cyan(genlog)
    fmt.Println(">")
    fmt.Println(">")
    fmt.Println(">")
    fmt.Println(">")


  }

  configdata := fmt.Sprintf("Initial population: %d Genes per organism:  %d Vertices per gene: %d Tournament size: %d Elite: %d Cutoff: %d Mutationrate: %f",population, polygons, vertices, tournamentsize, elite, cutoff, mutationrate)

  currentTime := time.Now() //currentTime.Format("2006-01-02-15:04:05")

  logpath = "logs/" + string(currentTime.Format("2006-01-02-15-04")) + "/"

  pathErr := os.MkdirAll(logpath ,0777)

	//check if you need to panic, fallback or report
	if pathErr != nil {
		fmt.Println(pathErr)
	}

  path := logpath + "config.txt"

  f, err := os.Create(path)
   if err != nil {
       fmt.Println(err)
       return
   }
   l, err := f.WriteString(configdata)
   if err != nil {
       fmt.Println(err)
       f.Close()
       return
   }
   fmt.Println(l, "configdata written successfully")
   err = f.Close()
   if err != nil {
       fmt.Println(err)
       return
   }




}

func waitStart(){

  buf := bufio.NewReader(os.Stdin)
   fmt.Print("> ")
   sentence, err := buf.ReadBytes('\n')
   if err != nil {
       fmt.Println(err)
   } else {
       fmt.Println(string(sentence))
   }

}

func stringint (str string) int {

  i1, err := strconv.Atoi(str)
  if err == nil {
       return i1
  }

  return i1
}


func loadReferenceImage(){

        imgfile, err := os.Open("ref.png")

        if err != nil {
                fmt.Println("img.jpg file not found!")
                os.Exit(1)
        }

        defer imgfile.Close()

        imgCfg, _, err := image.DecodeConfig(imgfile)

        if err != nil {
                fmt.Println(err)
                os.Exit(1)
        }

        width := imgCfg.Width
        height := imgCfg.Height

        xsize = width
        ysize = height

        fmt.Println(">")
        fmt.Println(">")
        fmt.Println(">")
        fmt.Println("Width : ", width)
        fmt.Println("Height : ", height)


        imgfile.Seek(0, 0)

        img, _, err := image.Decode(imgfile)

        var pixelsint []int

        for y := height - 1; y > -1; y-- {

                for x := 0; x < width; x++ {
                        r, g, b, a := img.At(x, y).RGBA()
                        pixelsint = append(pixelsint, int(float64(int(r)) / 65535 * 255))
                        pixelsint = append(pixelsint, int(float64(int(g)) / 65535 * 255))
                        pixelsint = append(pixelsint, int(float64(int(b)) / 65535 * 255))
                        pixelsint = append(pixelsint, int(float64(int(a)) / 65535 * 255))
                }

        }
        refimage = pixelsint


}


func loadPicture(path string) (pixel.Picture, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	img, _, err := image.Decode(file)
	if err != nil {
		return nil, err
	}
	return pixel.PictureDataFromImage(img), nil
}

func drawChart() {

  xvalues := make([]float64, Generation)
  for i:=0;i < Generation; i++ {

    xvalues[i] = float64(i)

  }

	graph := chart.Chart{
		XAxis: chart.XAxis{
			Name:      "Generation",
			NameStyle: chart.StyleShow(),
			Style:     chart.StyleShow(),
		},
		YAxis: chart.YAxis{
			Name:      "Fitness",
			NameStyle: chart.StyleShow(),
			Style:     chart.StyleShow(),
		},
		Series: []chart.Series{
			chart.ContinuousSeries{
				Style: chart.Style{
					Show:        true,
          DotColor: chart.GetDefaultColor(0).WithAlpha(64),
					StrokeColor: chart.GetDefaultColor(0).WithAlpha(64),
					FillColor:   chart.GetDefaultColor(0).WithAlpha(64),
				},
				XValues: xvalues,
				YValues: Fitnesses,
			},
		},
	}

  f, err := os.Create("chart.png")
  if err != nil {
    panic(err)
  }
  defer f.Close()

	graph.Render(chart.PNG, f)

}

func main() {
  rand.Seed(time.Now().UTC().UnixNano())
	pixelgl.Run(run)
}

var fs = require("fs")
var path = require("path")
var webshot = require("webshot")
var args = require('args')
var connect = require('connect')
var serveStatic = require('serve-static')
var getPort = require("get-port")

// Parse command line arguments
args
  .option("file", "Which JSON file to load.")
  .option("output", "Where to putput the picture.")

var flags = args.parse(process.argv)

// Check if file path is passed
if (!flags.file) {
  args.showHelp()
  process.exit(1)
}
// Read file as JSON and pass through to PhantomJS session for the snapshot.
getPort({ port: 2337 }).then(function(port) {
  fs.readFile(flags.file, 'utf8', function(err, data) {
    if (err) {
      console.error("Error reading json from: " + flags.file)
      process.exit(1)
    }
    var jsonFile = JSON.parse(data)
    var options = {
      renderDelay: 200,
      streamType: "jpg",
      captureSelector: "body",
      errorIfJSException: true,
      onInitialized: {
        // This function is invoked in the phantomJS session
        fn: function() { window.jsonFile = this.jsonFile },
        context: { jsonFile: jsonFile }
      }
    }
    // Take a snapshot of the plotting page running in PhantomJS.
    function chdir(dir) {
      try {
        process.chdir(dir); // go to webpage directory
      } catch (err) {
        console.error(err)
        process.exit(1)
      }
    }
    var oldCwd = process.cwd()
    var output = path.join(oldCwd, flags.output || "dependency-tree-image.jpg")
    chdir(__dirname)
    console.log(process.cwd())
    connect().use(serveStatic(__dirname)).listen(port, function() {
      console.log('Snapshot server running on port: ' + port);
      var url = "http://localhost:" + port + "/circular.html";
      webshot(url, output, options, function (error) {
        chdir(oldCwd); // go back to working directory
        if (error) {
          console.error(error);
          process.exit(1);
        }
        process.exit(0);
      })
    })
  })
}) 


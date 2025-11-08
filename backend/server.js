import express from "express";
import { exec } from "child_process";
import fs from "fs";
import path from "path";
import cors from "cors";

const app = express();
const PORT = 5000;

app.use(cors());
app.use(express.json());

app.get("/run-stress-test", (req, res) => {
  console.log("Running stress test...");

  // Execute compiled C++ binary (make sure it’s compiled first)
  exec("./stress_test", (error, stdout, stderr) => {
    if (error) {
      console.error(`Error executing stress test: ${error.message}`);
      return res.status(500).json({ error: "Error running stress test" });
    }

    // Read results.csv after program finishes
    const filePath = path.join(process.cwd(), "results.csv");
    fs.readFile(filePath, "utf8", (err, data) => {
      if (err) {
        console.error("Error reading results.csv", err);
        return res.status(500).json({ error: "Cannot read results.csv" });
      }

      // Parse CSV → JSON
      const lines = data.trim().split("\n");
      const headers = lines.shift().split(",");
      const json = lines.map(line => {
        const values = line.split(",");
        const obj = {};
        headers.forEach((h, i) => obj[h] = values[i]);
        return obj;
      });

      console.log("Stress test completed successfully.");
      res.json({ output: stdout, results: json });
    });
  });
});

app.listen(PORT, () => console.log(`✅ Server running on port ${PORT}`));

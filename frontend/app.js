const { useState } = React;

function App() {
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState([]);

  const runStressTest = async () => {
    setLoading(true);
    try {
      const res = await axios.get("http://localhost:5000/run-stress-test");
      setResults(res.data.results);
    } catch (err) {
      alert("Error running stress test");
      console.error(err);
    }
    setLoading(false);
  };

  return (
    <div>
      <h1>📊 Portfolio Stress Test Dashboard</h1>
      <button onClick={runStressTest} disabled={loading}>
        {loading ? "Running..." : "Run Stress Test"}
      </button>

      {results.length > 0 && (
        <div>
          <h2>Results</h2>
          <table>
            <thead>
              <tr>
                <th>Scenario</th>
                <th>Shock (%)</th>
                <th>Final Value</th>
              </tr>
            </thead>
            <tbody>
              {results.map((row, i) => (
                <tr key={i}>
                  <td>{row.Scenario}</td>
                  <td>{row["Shock(%)"]}</td>
                  <td>{parseFloat(row.FinalValue).toFixed(2)}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}

// Mount React app
const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(<App />);

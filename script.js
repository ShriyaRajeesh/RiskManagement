document.getElementById("fileInput").addEventListener("change", function (e) {
  const file = e.target.files[0];
  if (!file) return;

  const reader = new FileReader();
  reader.onload = function (event) {
    const text = event.target.result;
    const rows = text.trim().split("\n").slice(1);
    const scenarios = [], shocks = [], values = [];

    rows.forEach(row => {
      const [scenario, shock, value] = row.split(",");
      scenarios.push(`Scenario ${scenario}`);
      shocks.push(parseFloat(shock));
      values.push(parseFloat(value));
    });

    renderChart(scenarios, shocks, values);
  };
  reader.readAsText(file);
});

function renderChart(scenarios, shocks, values) {
  const ctx = document.getElementById("stressChart").getContext("2d");

  if (window.stressChart) {
    window.stressChart.destroy();
  }

  window.stressChart = new Chart(ctx, {
    type: "bar",
    data: {
      labels: scenarios,
      datasets: [{
        label: "Final Portfolio Value",
        data: values,
        backgroundColor: shocks.map(s => s < 0 ? "#DA0037" : "#42f5e6"),
      }]
    },
    options: {
      plugins: {
        title: {
          display: true,
          text: "Portfolio Value under Different Market Shocks",
          color: "#EDEDED",
          font: { size: 18 }
        },
        legend: { labels: { color: "#EDEDED" } }
      },
      scales: {
        x: { ticks: { color: "#EDEDED" } },
        y: {
          ticks: { color: "#EDEDED" },
          title: { display: true, text: "Final Portfolio Value ($)", color: "#EDEDED" }
        }
      }
    }
  });

  // Add summary info
  const minValue = Math.min(...values);
  const maxValue = Math.max(...values);
  document.getElementById("summary").innerText =
    `Lowest Value: $${minValue.toFixed(2)} | Highest Value: $${maxValue.toFixed(2)}`;
}

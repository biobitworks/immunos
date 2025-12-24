---
source: /Users/byron/projects/rockbeatspaper/src/routes/analysis/index.tsx
relative: rockbeatspaper/src/routes/analysis/index.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { createFileRoute } from "@tanstack/react-router";
import { useMutation, useQuery } from "@tanstack/react-query";
import { useState } from "react";
import { useForm } from "react-hook-form";
import toast from "react-hot-toast";
import { Brain, Sparkles, TrendingDown, AlertCircle } from "lucide-react";
import { ProtectedLayout } from "~/components/ProtectedLayout";
import { Card } from "~/components/Card";
import { Button } from "~/components/Button";
import { Input } from "~/components/Input";
import { useTRPC } from "~/trpc/react";
import { useAuthStore } from "~/stores/auth";

export const Route = createFileRoute("/analysis/")({
  component: Analysis,
});

type AnalysisResult = {
  analysisId: number;
  summary: string;
  recommendations: Array<{
    title: string;
    description: string;
    priority: "high" | "medium" | "low";
  }>;
  estimatedSavings: number;
  insights: string[];
};

function Analysis() {
  const trpc = useTRPC();
  const { authToken } = useAuthStore();
  const [analysisResult, setAnalysisResult] = useState<AnalysisResult | null>(
    null
  );

  const {
    register,
    handleSubmit,
    formState: { errors },
  } = useForm<{
    currentTemperature: number;
    targetTemperature: number;
    humidity: number;
    deviceType: string;
    location: string;
    currentPowerUsage: number;
  }>();

  const historyQuery = useQuery({
    ...trpc.getAnalysisHistory.queryOptions({ authToken: authToken! }),
    enabled: !!authToken,
  });

  const analyzeMutation = useMutation(
    trpc.analyzeThermalData.mutationOptions({
      onSuccess: (data) => {
        setAnalysisResult(data);
        toast.success("Analysis complete!");
      },
      onError: (error) => {
        toast.error(error.message || "Analysis failed");
      },
    })
  );

  const onSubmit = (data: {
    currentTemperature: number;
    targetTemperature: number;
    humidity: number;
    deviceType: string;
    location: string;
    currentPowerUsage: number;
  }) => {
    analyzeMutation.mutate({
      authToken: authToken!,
      currentTemperature: Number(data.currentTemperature),
      targetTemperature: Number(data.targetTemperature),
      humidity: Number(data.humidity),
      deviceType: data.deviceType,
      location: data.location,
      currentPowerUsage: Number(data.currentPowerUsage),
    });
  };

  const priorityColors = {
    high: "bg-red-100 text-red-800 border-red-200",
    medium: "bg-yellow-100 text-yellow-800 border-yellow-200",
    low: "bg-green-100 text-green-800 border-green-200",
  };

  return (
    <ProtectedLayout>
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="mb-8">
          <h1 className="text-3xl font-bold text-gray-900 dark:text-white mb-2">
            AI Thermal Analysis
          </h1>
          <p className="text-gray-600 dark:text-gray-400">
            Get intelligent recommendations to optimize your energy efficiency
          </p>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Analysis Form */}
          <Card>
            <div className="flex items-center mb-6">
              <div className="w-12 h-12 bg-gradient-to-r from-blue-600 to-purple-600 rounded-lg flex items-center justify-center mr-4">
                <Brain className="w-6 h-6 text-white" />
              </div>
              <div>
                <h2 className="text-xl font-bold text-gray-900 dark:text-white">
                  Input Data
                </h2>
                <p className="text-sm text-gray-600 dark:text-gray-400">
                  Enter your current thermal readings
                </p>
              </div>
            </div>

            <form onSubmit={handleSubmit(onSubmit)} className="space-y-4">
              <div className="grid grid-cols-2 gap-4">
                <Input
                  label="Current Temp (°C)"
                  type="number"
                  step="0.1"
                  {...register("currentTemperature", {
                    required: "Required",
                    valueAsNumber: true,
                  })}
                  error={errors.currentTemperature?.message}
                />
                <Input
                  label="Target Temp (°C)"
                  type="number"
                  step="0.1"
                  {...register("targetTemperature", {
                    required: "Required",
                    valueAsNumber: true,
                  })}
                  error={errors.targetTemperature?.message}
                />
              </div>

              <Input
                label="Humidity (%)"
                type="number"
                step="1"
                {...register("humidity", {
                  required: "Required",
                  valueAsNumber: true,
                })}
                error={errors.humidity?.message}
              />

              <Input
                label="Device Type"
                {...register("deviceType", { required: "Required" })}
                error={errors.deviceType?.message}
                placeholder="e.g., HVAC, Thermostat"
              />

              <Input
                label="Location"
                {...register("location", { required: "Required" })}
                error={errors.location?.message}
                placeholder="e.g., Living Room"
              />

              <Input
                label="Power Usage (W)"
                type="number"
                step="1"
                {...register("currentPowerUsage", {
                  required: "Required",
                  valueAsNumber: true,
                })}
                error={errors.currentPowerUsage?.message}
              />

              <Button
                type="submit"
                disabled={analyzeMutation.isPending}
                className="w-full"
              >
                {analyzeMutation.isPending ? (
                  <>Analyzing...</>
                ) : (
                  <>
                    <Sparkles className="w-4 h-4 mr-2" />
                    Analyze with AI
                  </>
                )}
              </Button>
            </form>
          </Card>

          {/* Analysis Results */}
          <div className="space-y-6">
            {analysisResult ? (
              <>
                <Card>
                  <div className="flex items-center justify-between mb-4">
                    <h2 className="text-xl font-bold text-gray-900 dark:text-white">
                      Analysis Results
                    </h2>
                    <div className="flex items-center space-x-2 text-green-600">
                      <TrendingDown className="w-5 h-5" />
                      <span className="text-2xl font-bold">
                        {analysisResult.estimatedSavings}%
                      </span>
                    </div>
                  </div>
                  <p className="text-gray-700 dark:text-gray-300 mb-4">
                    {analysisResult.summary}
                  </p>
                  <div className="bg-blue-50 dark:bg-blue-900/20 rounded-lg p-4">
                    <p className="text-sm font-semibold text-blue-900 dark:text-blue-200 mb-2">
                      Key Insights:
                    </p>
                    <ul className="space-y-1">
                      {analysisResult.insights.map((insight, idx) => (
                        <li
                          key={idx}
                          className="text-sm text-blue-800 dark:text-blue-300"
                        >
                          • {insight}
                        </li>
                      ))}
                    </ul>
                  </div>
                </Card>

                <Card>
                  <h3 className="text-lg font-bold text-gray-900 dark:text-white mb-4">
                    Recommendations
                  </h3>
                  <div className="space-y-3">
                    {analysisResult.recommendations.map((rec, idx) => (
                      <div
                        key={idx}
                        className={`border rounded-lg p-4 ${priorityColors[rec.priority]}`}
                      >
                        <div className="flex items-start justify-between mb-2">
                          <h4 className="font-semibold">{rec.title}</h4>
                          <span className="text-xs uppercase font-bold px-2 py-1 rounded">
                            {rec.priority}
                          </span>
                        </div>
                        <p className="text-sm">{rec.description}</p>
                      </div>
                    ))}
                  </div>
                </Card>
              </>
            ) : (
              <Card>
                <div className="text-center py-12">
                  <AlertCircle className="w-12 h-12 text-gray-400 mx-auto mb-4" />
                  <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                    No Analysis Yet
                  </h3>
                  <p className="text-gray-600 dark:text-gray-400">
                    Fill out the form to get AI-powered recommendations
                  </p>
                </div>
              </Card>
            )}
          </div>
        </div>

        {/* Analysis History */}
        {historyQuery.data?.analyses && historyQuery.data.analyses.length > 0 && (
          <Card className="mt-8">
            <h2 className="text-xl font-bold text-gray-900 dark:text-white mb-4">
              Recent Analyses
            </h2>
            <div className="space-y-3">
              {historyQuery.data.analyses.map((analysis) => (
                <div
                  key={analysis.id}
                  className="border border-gray-200 dark:border-gray-700 rounded-lg p-4 hover:shadow-md transition-shadow"
                >
                  <div className="flex items-center justify-between">
                    <div>
                      <p className="font-semibold text-gray-900 dark:text-white">
                        {analysis.analysisType}
                      </p>
                      <p className="text-sm text-gray-600 dark:text-gray-400">
                        {new Date(analysis.createdAt).toLocaleString()}
                      </p>
                    </div>
                    {analysis.potentialSavings && (
                      <div className="text-right">
                        <p className="text-sm text-gray-600 dark:text-gray-400">
                          Potential Savings
                        </p>
                        <p className="text-2xl font-bold text-green-600">
                          {analysis.potentialSavings}%
                        </p>
                      </div>
                    )}
                  </div>
                </div>
              ))}
            </div>
          </Card>
        )}
      </div>
    </ProtectedLayout>
  );
}

```

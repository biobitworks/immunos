---
source: /Users/byron/projects/rockbeatspaper/src/routes/pi-benford/index.tsx
relative: rockbeatspaper/src/routes/pi-benford/index.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { createFileRoute } from "@tanstack/react-router";
import { useMutation } from "@tanstack/react-query";
import { useState } from "react";
import { useForm } from "react-hook-form";
import toast from "react-hot-toast";
import { Calculator, Pi, TrendingUp, AlertCircle } from "lucide-react";
import { Card } from "~/components/Card";
import { Button } from "~/components/Button";
import { Input } from "~/components/Input";
import { useTRPC } from "~/trpc/react";

export const Route = createFileRoute("/pi-benford/")({
  component: PiBenfordAnalysis,
});

function PiBenfordAnalysis() {
  const trpc = useTRPC();
  const [analysisOutput, setAnalysisOutput] = useState<string | null>(null);

  const {
    register,
    handleSubmit,
    formState: { errors },
  } = useForm<{
    numDigits: number;
  }>({
    defaultValues: {
      numDigits: 1000,
    },
  });

  const analysisMutation = useMutation(
    trpc.runPiBenfordAnalysis.mutationOptions({
      onSuccess: (data) => {
        setAnalysisOutput(data.output);
        toast.success("Analysis complete!");
      },
      onError: (error) => {
        toast.error(error.message || "Analysis failed");
      },
    })
  );

  const onSubmit = (data: { numDigits: number }) => {
    analysisMutation.mutate({
      numDigits: Number(data.numDigits),
    });
  };

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-purple-50 dark:from-gray-900 dark:to-gray-800">
      <div className="max-w-6xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        {/* Header */}
        <div className="mb-8 text-center">
          <div className="flex items-center justify-center mb-4">
            <div className="w-16 h-16 bg-gradient-to-r from-blue-600 to-purple-600 rounded-full flex items-center justify-center">
              <Pi className="w-8 h-8 text-white" />
            </div>
          </div>
          <h1 className="text-4xl font-bold text-gray-900 dark:text-white mb-2">
            Pi Benford's Law Analysis
          </h1>
          <p className="text-gray-600 dark:text-gray-400 max-w-2xl mx-auto">
            Analyze the leading digits of pi's decimal expansion using Benford's Law
            and THRML (Thermodynamic Machine Learning)
          </p>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Input Form */}
          <Card>
            <div className="flex items-center mb-6">
              <div className="w-12 h-12 bg-gradient-to-r from-blue-600 to-purple-600 rounded-lg flex items-center justify-center mr-4">
                <Calculator className="w-6 h-6 text-white" />
              </div>
              <div>
                <h2 className="text-xl font-bold text-gray-900 dark:text-white">
                  Analysis Parameters
                </h2>
                <p className="text-sm text-gray-600 dark:text-gray-400">
                  Configure your analysis
                </p>
              </div>
            </div>

            <form onSubmit={handleSubmit(onSubmit)} className="space-y-6">
              <div>
                <Input
                  label="Number of Digits"
                  type="number"
                  {...register("numDigits", {
                    required: "Number of digits is required",
                    valueAsNumber: true,
                    min: {
                      value: 10,
                      message: "Minimum 10 digits required",
                    },
                    max: {
                      value: 100000,
                      message: "Maximum 100,000 digits allowed",
                    },
                  })}
                  error={errors.numDigits?.message}
                  placeholder="e.g., 1000"
                />
                <p className="mt-2 text-sm text-gray-600 dark:text-gray-400">
                  Enter the number of digits after the decimal point for π (between 10 and 100,000)
                </p>
              </div>

              <div className="bg-blue-50 dark:bg-blue-900/20 rounded-lg p-4">
                <h3 className="text-sm font-semibold text-blue-900 dark:text-blue-200 mb-2">
                  What is Benford's Law?
                </h3>
                <p className="text-sm text-blue-800 dark:text-blue-300">
                  Benford's Law states that in many naturally occurring datasets, 
                  the leading digit is more likely to be small. For example, 
                  the digit 1 appears as the leading digit about 30% of the time.
                </p>
              </div>

              <Button
                type="submit"
                disabled={analysisMutation.isPending}
                className="w-full"
              >
                {analysisMutation.isPending ? (
                  <>Running Analysis...</>
                ) : (
                  <>
                    <TrendingUp className="w-4 h-4 mr-2" />
                    Run Analysis
                  </>
                )}
              </Button>
            </form>
          </Card>

          {/* Results */}
          <div className="space-y-6">
            {analysisOutput ? (
              <Card>
                <h2 className="text-xl font-bold text-gray-900 dark:text-white mb-4">
                  Analysis Results
                </h2>
                <div className="bg-gray-900 dark:bg-gray-950 rounded-lg p-4 overflow-auto">
                  <pre className="text-sm text-green-400 font-mono whitespace-pre-wrap">
                    {analysisOutput}
                  </pre>
                </div>
              </Card>
            ) : (
              <Card>
                <div className="text-center py-12">
                  <AlertCircle className="w-12 h-12 text-gray-400 mx-auto mb-4" />
                  <h3 className="text-lg font-semibold text-gray-900 dark:text-white mb-2">
                    No Analysis Yet
                  </h3>
                  <p className="text-gray-600 dark:text-gray-400">
                    Enter the number of digits and run the analysis to see results
                  </p>
                </div>
              </Card>
            )}

            {/* Info Card */}
            <Card>
              <h3 className="text-lg font-bold text-gray-900 dark:text-white mb-3">
                About This Analysis
              </h3>
              <div className="space-y-3 text-sm text-gray-700 dark:text-gray-300">
                <p>
                  This analysis extracts leading digits from π's decimal expansion
                  and tests whether they follow Benford's Law using:
                </p>
                <ul className="list-disc list-inside space-y-1 ml-2">
                  <li>Chi-square goodness-of-fit test</li>
                  <li>THRML (Thermodynamic Machine Learning) categorical EBM</li>
                  <li>Gibbs sampling for frequency estimation</li>
                </ul>
                <p className="pt-2 border-t border-gray-200 dark:border-gray-700">
                  The analysis compares observed frequencies against theoretical
                  Benford's Law probabilities to determine if π's digits follow
                  the expected distribution.
                </p>
              </div>
            </Card>
          </div>
        </div>
      </div>
    </div>
  );
}

```

module.exports = {
		  branchPrefix: 'renovate/',
		  dryRun: false,
		  username: 'renovate-release',
		  gitAuthor: 'Renovate Bot <bot@renovateapp.com>',
		  onboarding: false,
		  platform: 'github',
		  requireConfig: false,
	          repositories: ["metagenomics/metagenomics-tk"],
		prHourlyLimit: 10,
  		prConcurrentLimit: 50,
  		packageRules: [
    {
      matchDatasources: ['docker'],
      matchPackagePrefixes: ['quay.io/biocontainers/'],
      versioning: 'regex:^(?<major>\\d+)\\.(?<minor>\\d+)\\.(?<patch>\\d+)--[^_]+_(?<build>\\d+)$'
    },
    {
      matchDatasources: ['docker'],
      excludePackagePrefixes: ['quay.io/biocontainers/'],
      versioning: 'loose'
    }
  ],
  regexManagers: [
    {
      fileMatch: ['^nextflow\\.config$'],
      matchStrings: [
        '[A-Za-z0-9_]+_image\\s*=\\s*configureImagePrefix\\("(?<depName>[^":]+(?:/[^":]+)+):(?<currentValue>[^"]+)"\\)'
      ],
      datasourceTemplate: 'docker'
    }
  ]
};
